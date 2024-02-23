#! /bin/sh
VERSION=latest

#
# Build local docker images
#
# When using docker the pipleine can use local images or pull them from DockerHub. 
# If you want to build the images yourself you can do it like this:
#
#   cd envs
#   ./build-docker-imgs.sh
#

docker build --platform=linux/amd64 -t haniffalab/nf-soupx:${VERSION} -f ./soupx/Dockerfile .
docker build --platform=linux/amd64 -t haniffalab/nf-scanpy:${VERSION} -f ./scanpy/Dockerfile .
