// Global variables
#define BLOCK_SIZE_1D 512
#define BLOCK_SIZE_2D 16
#define WARP_SIZE 32
#define TILE_SIZE 32
#define WPT_X 4
#define WPT_Y 4
#define SHARED_MEMORY 49152
#ifdef USE_DOUBLES
typedef double prec_typ;
#else
typedef float prec_typ;
#endif
