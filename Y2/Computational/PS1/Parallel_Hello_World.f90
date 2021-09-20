PROGRAM Parallel_Hello_World
USE OMP_LIB

PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()

END