#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "defs.h"
#include <json-c/json.h>


// Write a single node as JSON object
void write_node(FILE *fp, const ColumnChunkMeta *node) {
    fprintf(fp, "{\n");
    fprintf(fp, "  \"rec_num\": %llu,\n", (unsigned long long)node->rec_num);
    fprintf(fp, "  \"file_offset\": %llu,\n", (unsigned long long)node->file_offset);
    fprintf(fp, "  \"uncompressed_size\": %llu,\n", (unsigned long long)node->uncompressed_size);
    fprintf(fp, "  \"compressed_size\": %llu\n", (unsigned long long)node->compressed_size);
    fprintf(fp, "}");
}

// Write a linked list as JSON array
void write_list(FILE *fp, ColumnChunkMeta *head) {
    fprintf(fp, "[\n");
    ColumnChunkMeta *current = head;
    while(current->prev) {
        current = current->prev; // Move to the head of the list
    }
    while (current)
    {
        write_node(fp, current);
        current = current->next;
        fprintf(fp, current ? ",\n" : "\n");
    }
    fprintf(fp, "]");
}

// Write entire array of linked lists as JSON object
void write_meta(FILE *fp, ColumnChunkMeta **array, size_t size) {
    fprintf(fp, "{\n");
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "  \"%s\": ", ColumnTypeNames[i]);
        write_list(fp, array[i]);
        fprintf(fp, "%s\n", i < size - 1 ? "," : "");
    }
    fprintf(fp, "}\n");
}

void write_header(FILE *fp, int64_t seekpos, int64_t meta_size) {
    fprintf(fp, "{");
    fprintf(fp, "\"seekpos\": %lld,\n", (long long)seekpos);
    fprintf(fp, "\"meta_size\": %lld\n", (long long)meta_size);
    fprintf(fp, "}\0");
}

void parse_meta_from_json_string(char *json_str, Reader *reader) {
    ColumnChunkMeta **array = (ColumnChunkMeta **)malloc(COLUMNTYPE_SIZE * sizeof(ColumnChunkMeta *));
    
    struct json_object * root = json_tokener_parse(json_str);

    if (root == NULL) {
        fprintf(stderr, "Failed to parse JSON string\n");
        exit(1);
    }

    for (size_t i = 0; i < COLUMNTYPE_SIZE; i++)
    {
        array[i] = NULL; // Initialize each pointer to NULL
        struct json_object *field = NULL; 
        if (!json_object_object_get_ex(root, ColumnTypeNames[i], &field)) {
            fprintf(stderr, "Key '%s' not found in JSON\n", ColumnTypeNames[i]);
            continue;
        }
        int64_t arr_size = json_object_array_length(field);
        array[i] = (ColumnChunkMeta *)malloc(arr_size * sizeof(ColumnChunkMeta));
        for (int64_t j = 0; j < arr_size; j++) {
            struct json_object *node = json_object_array_get_idx(field, j);
            struct json_object *rec_num, *file_offset_obj, *uncompressed_size_obj, *compressed_size_obj;

            json_object_object_get_ex(node, "rec_num", &rec_num);
            json_object_object_get_ex(node, "file_offset", &file_offset_obj);
            json_object_object_get_ex(node, "uncompressed_size", &uncompressed_size_obj);
            json_object_object_get_ex(node, "compressed_size", &compressed_size_obj);

            array[i][j].rec_num = json_object_get_int64(rec_num);
            array[i][j].file_offset = json_object_get_int64(file_offset_obj);
            array[i][j].uncompressed_size = json_object_get_int64(uncompressed_size_obj);
            array[i][j].compressed_size = json_object_get_int64(compressed_size_obj);
            array[i][j].next = NULL;
            array[i][j].prev = NULL;
        }
        reader->metadatas_lengths[i] = arr_size;
    }

    json_object_put(root); // Free the JSON object

    reader->metadatas = array;

    return array;
}