


# Script usage.
# build_macos [source_dir] [build_type] [build_directory] [install_directory]


source_dir=""
build_type=""  # Default to empty package
build_dir=""  # Default to empty target
build_install=""

# Parse options
while getopts ":s:t:b:p:" opt; do
  case ${opt} in
    h )
      echo "Usage:"
      echo "    pip -h                      Display this help message."
      echo "    pip install <package>       Install <package>."
      exit 0
      ;;
   \? )
     echo "Invalid Option: -$OPTARG" 1>&2
     exit 1
     ;;
  esac
done
shift $((OPTIND -1))

subcommand=$1; shift  # Remove 'pip' from the argument list
case "$subcommand" in
  # Parse options to the install sub command
  install)
    package=$1; shift  # Remove 'install' from the argument list

    # Process package options
    while getopts ":t:" opt; do
      case ${opt} in
        t )
          target=$OPTARG
          ;;
        \? )
          echo "Invalid Option: -$OPTARG" 1>&2
          exit 1
          ;;
        : )
          echo "Invalid Option: -$OPTARG requires an argument" 1>&2
          exit 1
          ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
esac

