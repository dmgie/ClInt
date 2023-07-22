def docker_current_dir() {
    REF_FULLPATH = "realpath ${params.reference_file}".execute().text.trim()
    "docker run -v $PWD:$PWD -v \$PWD:\$PWD -v ${REF_FULLPATH}:${REF_FULLPATH} --user ${params.USER_ID}:${params.GROUP_ID}"
}
