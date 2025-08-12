#ifndef GNB_FRAME_SLOT_BUFFER_H
#define GNB_FRAME_SLOT_BUFFER_H

#define NUM_FRAMES 1024
#define THRESHOLD 5
#define MAX_BUF_SLOTS 10 

#include <stdbool.h>
#include <pthread.h>
#include <stdint.h>

typedef struct {
  int frame;
  int slot;
} frame_slot;

extern frame_slot SLOT_BUFFER[THRESHOLD];
extern pthread_mutex_t SLOT_BUFFER_LOCK;

extern uint32_t slot_identifier;
#endif