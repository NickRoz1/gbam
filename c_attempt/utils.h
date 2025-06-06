#include <stdint.h>

// Write 16-bit integer to buffer (little-endian)
void write_int16_le(uint8_t *buf, int16_t value) {
    buf[0] = value & 0xFF;         // LSB first
    buf[1] = (value >> 8) & 0xFF;  // Next byte
}

// Write 32-bit integer to buffer (little-endian)
void write_int32_le(uint8_t *buf, int32_t value) {
    buf[0] = value & 0xFF;          // Byte 0 (LSB)
    buf[1] = (value >> 8) & 0xFF;   // Byte 1
    buf[2] = (value >> 16) & 0xFF;  // Byte 2
    buf[3] = (value >> 24) & 0xFF;  // Byte 3 (MSB)
}

// Read 16-bit integer from buffer (little-endian)
int16_t read_int16_le(const uint8_t *buf) {
    return buf[0] | (buf[1] << 8);
}

// Read 32-bit integer from buffer (little-endian)
int32_t read_int32_le(const uint8_t *buf) {
    return buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
}
