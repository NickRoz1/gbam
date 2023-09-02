import json
from collections import OrderedDict
import lz4.block as lb
import struct
import sys
import functools 

f = open("test_data/little.gbam", "rb")
header = (f.read(1000).rstrip(b'\x00'))
header_json = json.loads(header.decode('utf-8'))
f.seek(header_json["seekpos"])
meta_json = json.loads(f.read().decode('utf-8'))

fields_mapping = OrderedDict()
fields_mapping["RefID"] = None
fields_mapping["Pos"] = None
fields_mapping["LName"] = None
fields_mapping["Mapq"] = None
fields_mapping["Bin"] = None
fields_mapping["NCigar"] = None
fields_mapping["Flags"] = None
fields_mapping["SequenceLength"] = None
fields_mapping["NextRefID"] = None
fields_mapping["NextPos"] = None
fields_mapping["TemplateLength"] = None
fields_mapping["RawSeqLen"] = None
fields_mapping["RawTagsLen"] = None
fields_mapping["ReadName"] = "LName"
fields_mapping["RawCigar"] = "NCigar"
fields_mapping["RawSequence"] = "RawSeqLen"
fields_mapping["RawQual"] = "SequenceLength"
fields_mapping["RawTags"] = "RawTagsLen"

@functools.cache
def get_decompressed(seekpos, block_size, uncompressed_size):
    f.seek(seekpos)
    compressed_buf = f.read(block_size)
    return lb.decompress(compressed_buf, uncompressed_size)

def get_block(i, meta):
    blocks = meta["blocks"]
    loaded = blocks[0]["numitems"]

    cur = 0
    while loaded <= i:
        cur += 1
        loaded += blocks[cur]["numitems"]
    
    dec =  get_decompressed(blocks[cur]["seekpos"], blocks[cur]["block_size"], blocks[cur]["uncompressed_size"])
    assert(len(dec) == blocks[cur]["uncompressed_size"])
    return (loaded-blocks[cur]["numitems"], dec)

class FixedColumn():
    def __init__(self, field):
        self.meta = meta_json["field_to_meta"][field]

    def get_item(self, i):
        loaded, uncompressed_buf = get_block(i, self.meta)

        i -= loaded
        offset = i*self.meta["item_size"]

        return uncompressed_buf[offset:(offset+self.meta["item_size"])]

class VarColumn():
    def __init__(self, fixed_index_column, field):
        self.meta = meta_json["field_to_meta"][field]
        self.index = fixed_index_column

    def get_item(self, i):
        loaded, uncompressed_buf = get_block(i, self.meta)
        prev_end = 0

        if i-loaded > 0:
            prev_end = struct.unpack('<L', self.index.get_item(i-1))[0]
        
        end = struct.unpack('<L', self.index.get_item(i))[0]
        
        return uncompressed_buf[prev_end:end]

records_num = 0
for block in meta_json["field_to_meta"]["Flags"]["blocks"]:
    records_num += block["numitems"]

columns = {}
for key, value in fields_mapping.items():
    if value is None:
        columns[key] = FixedColumn(key)
    else:
        index = FixedColumn(value)
        columns[key] = VarColumn(index, key)

sys.stdout.buffer.write(b"BAM\x01")
sys.stdout.buffer.write(bytes(meta_json["sam_header"]))

for rec_i in range(0, records_num):
    arr = []
    arr.append(columns["RefID"].get_item(rec_i))
    arr.append(columns["Pos"].get_item(rec_i))
    arr.append(len(columns["ReadName"].get_item(rec_i)).to_bytes(1, "little"))
    arr.append(columns["Mapq"].get_item(rec_i))
    arr.append(columns["Bin"].get_item(rec_i))
    arr.append(struct.pack('<H', len(columns["RawCigar"].get_item(rec_i))//4))
    arr.append(columns["Flags"].get_item(rec_i))
    arr.append(struct.pack('<L', len(columns["RawQual"].get_item(rec_i))))
    arr.append(columns["NextRefID"].get_item(rec_i))
    arr.append(columns["NextPos"].get_item(rec_i))
    arr.append(columns["TemplateLength"].get_item(rec_i))
    arr.append(columns["ReadName"].get_item(rec_i))
    arr.append(columns["RawCigar"].get_item(rec_i))
    # t = bytearray(columns["RawSequence"].get_item(rec_i))
    # if len(t) > 0 and t[-1]&0x0F == 0:
    #     t[-1] = t[-1]|15
    # arr.append(t)
    arr.append(columns["RawSequence"].get_item(rec_i))
    arr.append(columns["RawQual"].get_item(rec_i))
    arr.append(columns["RawTags"].get_item(rec_i))
    
    total_len = 0
    for a in arr:
        total_len += len(a)
    arr.insert(0, struct.pack('<L', total_len))

    for a in arr:
        sys.stdout.buffer.write(a)

# You can pipe this output into SAMTOOLS VIEW: python3 minimal_parser.py | samtools view