cmake_minimum_required(VERSION 3.20)
project(epcsaft_native_sdk_consumer LANGUAGES CXX)

find_package(epcsaft CONFIG REQUIRED)

add_executable(epcsaft_native_sdk_consumer main.cpp)
target_link_libraries(epcsaft_native_sdk_consumer PRIVATE epcsaft::native_sdk)
target_compile_features(epcsaft_native_sdk_consumer PRIVATE cxx_std_17)
