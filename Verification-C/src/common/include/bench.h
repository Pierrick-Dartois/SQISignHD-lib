#include <time.h>

static inline int64_t cpucycles(void){
	struct timespec time;
	clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
}