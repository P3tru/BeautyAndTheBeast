#!/bin/bash

echo "Setting up ROOT_INCLUDE_PATH"
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${PWD}/include

echo "Setting up LD_LIBRARY_PATH"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/lib
