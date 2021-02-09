# ---------------------------------------------------------------------------- #
# emDNA cli
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


add_subdirectory(app/)
add_subdirectory(dna_force_field_packager/)
add_subdirectory(forceprobe/)
add_subdirectory(optbinding/)
add_subdirectory(tools/)


# Global target to build all cli tools.
add_custom_target(build-emDNA-tools
                  DEPENDS
                  emDNA-cli
                  emDNA-ff-packager
                  emDNA-force-probe
                  emDNA-pro-bind
                  emDNA-check-collision
                  emDNA-lego-gen
                  emDNA-parser
                  emDNA-topology
                  COMMENT "Building all emDNA tools")