# Tool versions and Containers

> **Note**: This fork adds additional containers for CRISPR editing QC and somatic variant analysis. Core pipeline containers are inherited from the upstream PacBio pipeline.

Containers are used to package tools and their dependencies. This ensures that the tools are reproducible and can be run on any system that supports the container runtime. Our containers are built using [Docker](https://www.docker.com/) and are compatible with any container runtime that supports the OCI Image Specification, like [Singularity](https://sylabs.io/singularity/) or [Podman](https://podman.io/).

Most of our containers are built on the `pb_wdl_base` container, which includes common bioinformatics tools and libraries.  We tag our containers with a version number and build count, but the containers are referenced within the WDL files by the sha256 sum tags for reproducibility and better compatibility with Cromwell and miniwdl call caching.

Our Dockerfiles can be inspected on GitHub, and the containers can be pulled from our [Quay.io organization](https://quay.io/repository/pacbio).

We directly use `deepvariant`, `deepvariant-gpu`, `pharmcat`, and `glnexus` containers from their respective authors, although we have mirrored some for better compatibility with Cromwell call caching.

| Container | Major tool versions | Dockerfile | Container |
| --------: | ------------------- | :---: | :---: |
| pb_wdl_base | <ul><li>htslib 1.20</li><li>bcftools 1.20</li><li>samtools 1.20</li><li>bedtools 2.31.0</li><li>python3.9</li><li>numpy 1.24.24</li><li>pandas 2.0.3</li><li>matplotlib 3.7.5</li><li>seaborn 0.13.2</li><li>pysam 0.22.1</li><li>vcfpy 0.13.8</li><li>biopython 1.83</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/pb_wdl_base) | [sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87](https://quay.io/repository/pacbio/pb_wdl_base/manifest/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87) |
| pbmm2 | <ul><li>pbmm2 1.17.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/9591749da92ca57f7283ca1c2268789c45fa341d/docker/pbmm2) | [pbmm2@sha256:5f3f4d1f5dbea5cd4c388ee26b2fecbbb7dbcef449343633e039dca3d3725859](https://quay.io/repository/pacbio/pbmm2/manifest/sha256:5f3f4d1f5dbea5cd4c388ee26b2fecbbb7dbcef449343633e039dca3d3725859) |
| mosdepth | <ul><li>mosdepth 0.3.9</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/fa84fbf582738c05c750e667ff43d11552ad4183/docker/mosdepth) | [mosdepth@sha256:63f7a5d1a4a17b71e66d755d3301a951e50f6b63777d34dab3ee9e182fd7acb1](https://quay.io/repository/pacbio/mosdepth/manifest/sha256:63f7a5d1a4a17b71e66d755d3301a951e50f6b63777d34dab3ee9e182fd7acb1) |
| sawfish | <ul><li>sawfish 2.2.1</li><li>sawshark 0.3.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/a9e9414ca16b5b25443b4352603551871d5683f3/docker/sawfish) | [sawfish@sha256:18ba096219fea38d6b32f5706fb794a05cc5d1d6cc16e2a09e3a13d62d8181d4](https://quay.io/repository/pacbio/sawfish/manifest/sha256:18ba096219fea38d6b32f5706fb794a05cc5d1d6cc16e2a09e3a13d62d8181d4) |
| trgt | <ul><li>trgt 5.0.0</li><li>`/opt/scripts/find_trgt_dropouts.py` 0.3.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/d9c818f3547e8c33cc6a9f1a65e311ec26db8569/docker/trgt) | [trgt@sha256:be0ed7c173d221bd84e360b2b056e2abbecadd07ed86ffd4883a5cecca7a1e57](https://quay.io/repository/pacbio/trgt/manifest/sha256:be0ed7c173d221bd84e360b2b056e2abbecadd07ed86ffd4883a5cecca7a1e57) |
| hiphase | <ul><li>hiphase 1.5.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/69039c010ada793bab4d38a9bd17a30562b9b671/docker/hiphase) | [hiphase@sha256:353b4ffdae4281bdd5daf5a73ea3bb26ea742ef2c36e9980cb1f1ed524a07482](https://quay.io/repository/pacbio/hiphase/manifest/sha256:353b4ffdae4281bdd5daf5a73ea3bb26ea742ef2c36e9980cb1f1ed524a07482) |
| mitorsaw | <ul><li>mitorsaw 0.2.4</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6f6cf280c8ac0b76dd1d08bd830347b0b8ca9cea/docker/mitorsaw) | [mitorsaw@sha256:d0e47fb84e6e962f01a754d1052a24e550694646c0d4afb056c0e3fd7ace7a0d](https://quay.io/repository/pacbio/mitorsaw/manifest/sha256:d0e47fb84e6e962f01a754d1052a24e550694646c0d4afb056c0e3fd7ace7a0d) |
| paraphase | <ul><li>paraphase 3.4.0</li><li>minimap 2.28</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/f25a5f465af066496e955e642e284cc45f378a76/docker/paraphase) | [paraphase@sha256:7e70bbc6666a33af9253f2df15dbbd57a7a031d40b166a02b58bf003d9932c4c](https://quay.io/repository/pacbio/paraphase/manifest/sha256:7e70bbc6666a33af9253f2df15dbbd57a7a031d40b166a02b58bf003d9932c4c) |
| pbstarphase | <ul><li>pbstarphase 2.0.0</li><li>Database 20251106</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/c52ac389edb76ba9224e7e2710d202aa0a4b840d/docker/pbstarphase) | [pbstarphase@sha256:b63d7abb718a11c29080cef19026d7c38b8d87df957969dd17b27b15c56de68b](https://quay.io/repository/pacbio/pbstarphase/manifest/sha256:b63d7abb718a11c29080cef19026d7c38b8d87df957969dd17b27b15c56de68b) |
| pb-cpg-tools | <ul><li>pb-cpg-tools 3.0.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/330b99b79f32b2d2598e812779f3c64460739e6c/docker/pb-cpg-tools) | [pb-cpg-tools@sha256:afd5468a423fe089f1437d525fdc19c704296f723958739a6fe226caa01fba1c](https://quay.io/repository/pacbio/pb-cpg-tools/manifest/sha256:afd5468a423fe089f1437d525fdc19c704296f723958739a6fe226caa01fba1c) |
| methbat | <ul><li>methbat 0.15.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/72290493cb204e70be35479e250d83d5ec35df19/docker/methbat) | [methbat@sha256:54a74389cf8ac485e8f34522d48e05880c01245e9aaf4dec6a6eddf25ee2c550](https://quay.io/repository/pacbio/methbat/manifest/sha256:54a74389cf8ac485e8f34522d48e05880c01245e9aaf4dec6a6eddf25ee2c550) |
| wgs_tertiary | <ul><li>`/opt/scripts/calculate_phrank.py` 2.0.0</li><li>`/opt/scripts/json2ped.py` 0.5.0</li></ul>Last built 2021-09-17:<ul><li>ensembl -> HGNC</li><li>ensembl -> HPO</li><li>HGNC -> inheritance</li><li>HPO DAG</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/fd70e2872bd3c6bb705faff5bc68374116d7d62f/docker/wgs_tertiary) | [wgs_tertiary@sha256:410597030e0c85cf16eb27a877d260e7e2824747f5e8b05566a1aaa729d71136](https://quay.io/repository/pacbio/wgs_tertiary/manifest/sha256:410597030e0c85cf16eb27a877d260e7e2824747f5e8b05566a1aaa729d71136) |
| slivar | <ul><li>slivar 0.3.1</li><li>`/opt/scripts/add_comphet_phase.py` 0.1.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/5e1094fd6755203b4971fdac6dcb951bbc098bed/docker/slivar) | [slivar@sha256:f71a27f756e2d69ec30949cbea97c54abbafde757562a98ef965f21a28aa8eaa](https://quay.io/repository/pacbio/slivar/manifest/sha256:f71a27f756e2d69ec30949cbea97c54abbafde757562a98ef965f21a28aa8eaa) |
| svpack | <ul><li>svpack 54b54db</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6fc750b0c65b4a5c1eb65791eab9eed89864d858/docker/svpack) | [svpack@sha256:628e9851e425ed8044a907d33de04043d1ef02d4d2b2667cf2e9a389bb011eba](https://quay.io/repository/pacbio/svpack/manifest/sha256:628e9851e425ed8044a907d33de04043d1ef02d4d2b2667cf2e9a389bb011eba) |
| deepvariant | <ul><li>DeepVariant 1.9.0</li></ul> |  | [deepvariant:1.9.0](https://hub.docker.com/layers/google/deepvariant/1.9.0/images/sha256-07e95b34e40cc50074d23273d479934a27e80919ac75bd97bf39a731e3c2d6ad) |
| deepvariant-gpu | <ul><li>DeepVariant 1.9.0</li></ul> |  | [deepvariant:1.9.0-gpu](https://hub.docker.com/layers/google/deepvariant/1.9.0-gpu/images/sha256-e0c8734b8700d945e3ee78d609acb90548f829c874596ffca436af8cf379f87a) |
| pharmcat | <ul><li>PharmCat 2.15.4</li></ul> |  | [pharmcat:2.15.4](https://hub.docker.com/layers/pgkb/pharmcat/2.15.4/images/sha256-5b58ae959b4cd85986546c2d67e3596f33097dedc40dfe57dd845b6e78781eb6) |
| glnexus | <ul><li>GLnexus 1.4.3</li></ul> |  | [glnexus:1.4.3](https://quay.io/repository/pacbio/glnexus/manifest/sha256:ce6fecf59dddc6089a8100b31c29c1e6ed50a0cf123da9f2bc589ee4b0c69c8e) |

## Additional Containers for CRISPR Editing QC

This fork adds specialized tools for CRISPR editing validation and somatic variant analysis:

| Container | Major tool versions | Container |
| --------: | ------------------- | :---: |
| bcftools (biocontainers) | <ul><li>bcftools (alternate version)</li></ul> | [bcftools@sha256:065115d2fbdcba9b9a9d360fe9f024ba829de2e64309a637a274ed47e0af17db](https://quay.io/repository/biocontainers/bcftools) |
| truvari | <ul><li>truvari (parent filtering)</li></ul> | [truvari@sha256:022be7a2dfe5decf1aff8d007a914f56ed68d0bc224888e1764820922bc11935](https://quay.io/repository/biocontainers/truvari) |
| minimap2 | <ul><li>minimap2 (edit alignment)</li></ul> | [minimap2@sha256:fdc9ef8bfbd31bab59a61b2e90dd226647deed971556175a4cd004f0bcdc7608](https://quay.io/repository/biocontainers/minimap2) |
| guidescan | <ul><li>guidescan (off-target detection)</li></ul> | [guidescan@sha256:9d93021243780b1faff47f1df4c1ae495177ff65ccc8273f0ec590caad5c82f0](https://quay.io/repository/biocontainers/guidescan) |
| cnvpytor | <ul><li>CNVpytor (CNV plotting)</li></ul> | [cnvpytor@sha256:42e480d51b4c640ebb46ede4ded05486c9974c0d0bde3a3cef9884ad56674383](https://quay.io/repository/biocontainers/cnvpytor) |
| pybigwig | <ul><li>pyBigWig 0.3.24</li></ul> | [pybigwig:0.3.24](https://quay.io/repository/biocontainers/pybigwig) |

## Additional Containers for Somatic Variant Analysis

| Container | Major tool versions | Container |
| --------: | ------------------- | :---: |
| ensembl-vep | <ul><li>VEP (variant annotation)</li></ul> | [ensembl-vep@sha256:e7612ab7c2923f2b9a78592b939e74874cd29f7494d70ee7135c8303841b03a8](https://hub.docker.com/r/ensemblorg/ensembl-vep) |
| annotsv | <ul><li>AnnotSV (SV annotation)</li></ul> | [annotsv@sha256:0c73fef5fa529b11e10bea0355480f01b56d0feb21af54cb9bbbd1f9f4c862a7](https://quay.io/repository/biocontainers/annotsv) |
| severus | <ul><li>Severus (phased SV calling)</li></ul> | [severus@sha256:fb4471e0504d564de78215ae15c081a1bb2022ad51e993eba92bc6fa5052a05d](https://quay.io/repository/biocontainers/severus) |
| somatic_general_tools | <ul><li>Somatic analysis utilities</li></ul> | [somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6](https://quay.io/repository/pacbio/somatic_general_tools) |
| chord | <ul><li>CHORD (HRD prediction)</li></ul> | [chord@sha256:9f6aa44ffefe3f736e66a0e2d7941d4f3e1cc6d848a9a11a17e85a6525e63a77](https://hub.docker.com/r/scwatts/chord) |
