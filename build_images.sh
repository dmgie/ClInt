
docker build --tag clint/clint:latest - < Dockerfile
singularity build clint.sif docker-daemon://clint/clint:latest
