#!/usr/bin/env bash
set -e

# This script generates a manifest file that lists all the docker images used in the WDL files.

grep '@sha' -h -r workflows/ \
| tr --squeeze-repeats ' ' \
| cut --fields=3 --delimiter=' ' \
| sed 's!~{runtime_attributes.container_registry}!quay.io/pacbio!;s/"//g;' \
| sort --unique \
> ./image_manifest.txt

# Also capture tag-pinned images (docker: "repo:tag" with no @sha256 digest). These are
# otherwise missed by the @sha grep above, so they never get prepulled and pull live at
# runtime -- a common cause of build-temp/registry failures. Exclude digest-pinned images,
# the runtime_attributes registry interpolation, and the locally-built `knock-knock` image
# (built separately by scripts/setup.sh).
grep -hroE 'docker:[[:space:]]*"[^"]+"' workflows/ \
| sed -E 's/^docker:[[:space:]]*"//; s/"$//' \
| grep -vE '@sha256|~\{|^knock-knock$' \
| sort --unique \
>> ./image_manifest.txt

deepvariant_version=1.9.0
echo "google/deepvariant:${deepvariant_version}" >> ./image_manifest.txt
echo "google/deepvariant:${deepvariant_version}-gpu" >> ./image_manifest.txt

pharmcat_version=2.15.4
echo "pgkb/pharmcat:${pharmcat_version}" >> ./image_manifest.txt