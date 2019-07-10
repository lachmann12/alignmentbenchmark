#####################################
#   AUTHOR: ALEXANDER LACHMANN      #
#   DATE: 5/16/2019                 #
#   Mount Sinai School of Medicine  #
#####################################

now=$(date +"%T")
echo "Current time : $now"
echo "------------------------------------"
# get rid of old stuff
docker rmi -f $(docker images | grep "^<none>" | awk "{print $3}")
docker rm $(docker ps -q -f status=exited)
docker system -y prune

docker build --no-cache -f DockerAlign -t maayanlab/alignmentbenchmark .

docker push maayanlab/alignmentbenchmark
