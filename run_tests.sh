#!/bin/bash
# Build and run tests in Docker

set -e

echo "Building Docker image..."
docker build -t hex-tests .

echo ""
echo "Running tests..."
docker run --rm hex-tests

echo ""
echo "Tests completed successfully!"
