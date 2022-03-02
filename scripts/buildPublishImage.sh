COMMIT_START=$1
COMMIT_END=$2
DOCKER_REPOSITORY=$3

for versionFilePath in $(git diff-tree --no-commit-id --name-only -r ${COMMIT_END} ${COMMIT_START} | grep "VERSION");
do
  folder=${versionFilePath%"/VERSION"}
  IMAGE_NAME=${folder##*/}

  tmpName="image-$RANDOM"
  docker build $folder --file $folder/Dockerfile --tag $tmpName
  IMAGE_ID=${DOCKER_REPOSITORY}/$IMAGE_NAME
  VERSION=$(cat $versionFilePath)

  echo IMAGE_ID=$IMAGE_ID
  echo VERSION=$VERSION

  docker tag $tmpName $IMAGE_ID:$VERSION
  docker push $IMAGE_ID:$VERSION
done;

