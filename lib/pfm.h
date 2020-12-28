
/* Read/write 1 and 3 channel PFM files, public domain Connelly Barnes 2007. */

#ifndef _pfm_h
#define _pfm_h

float *read_pfm_file(const char *filename, unsigned *w, unsigned *h);
void write_pfm_file(const char *filename, float *depth, unsigned w, unsigned h);
float *read_pfm_file3(const char *filename, unsigned *w, unsigned *h);
void write_pfm_file3(const char *filename, float *depth, unsigned w, unsigned h);

#endif
