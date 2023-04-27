// Source: https://github.com/miloyip/itoa-benchmark

static G_DIGITS_LUT: [char; 200] = [
    '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
    '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
    '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
    '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
    '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
    '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
    '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
    '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
    '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
    '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9'
];

static POWERS_OF_10: [u32; 10] = [
    0,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000,
    1000000000
];

fn count_decimal_digit32(num: u32) -> usize {
    let t: usize = (((32 - (num | 1).leading_zeros()) * 1233) >> 12) as usize;
    t - (num < POWERS_OF_10[t]) as usize + 1
}
    
/// Additional count number of digit pass
/// Use lookup table of two G_DIGITS_LUT
/// # Safety
/// Make sure buffer can actually fit u32 integer
pub unsafe fn u32toa_countlut(mut value: u32, mut buffer: *mut u8) -> *mut u8 {
    let digit = count_decimal_digit32(value);
    unsafe {
        buffer = buffer.add(digit);

        // ptr is right after the number in the buffer, you can place space/tab here.
        while value >= 100 {
            let i = ((value % 100) << 1) as usize;
            value /= 100;
            buffer = buffer.offset(-1);
            *buffer = G_DIGITS_LUT[i + 1] as u8;
            buffer =  buffer.offset(-1);
            *buffer = G_DIGITS_LUT[i] as u8;
        }

        if value < 10 {
            buffer = buffer.offset(-1);
            *buffer = value as u8 + b'0';
        }
        else {
            let i = (value << 1) as usize;
            buffer = buffer.offset(-1);
            *buffer = G_DIGITS_LUT[i + 1] as u8;
            buffer = buffer.offset(-1);
            *buffer = G_DIGITS_LUT[i] as u8;
        }
        buffer.add(digit)
    }
}

/// # Safety
/// Make sure buffer can actually fit i32 integer
pub unsafe fn i32toa_countlut(value: i32, mut buffer: *mut u8) -> *mut u8 {
    let mut u:u32 = value as u32;
    if value < 0 {
        unsafe {
            *buffer = b'-';
            buffer = buffer.add(1);
        }
        u = !u + 1;
    }
    u32toa_countlut(u, buffer)
}