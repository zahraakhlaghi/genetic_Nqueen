cmake_minimum_required(VERSION 3.17)
project(ga C)

set(CMAKE_C_STANDARD 99)

add_executable(ga main.c eval.c gen.c init.c random.c rep.c select.c stats.c xover.c utils.c)

target_link_libraries(ga PRIVATE m)