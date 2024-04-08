include { INTEGRATE } from "../modules/integrate.nf"
include { INTEGRATE as INTEGRATE_GPU } from "../modules/integrate.nf"
include { INTEGRATE as INTEGRATE_SCVI } from "../modules/integrate.nf"
include { INTEGRATE_SCANVI } from "../modules/integrate_scanvi.nf"
include { INTEGRATE_SCARCHES } from "../modules/integrate_scarches.nf"
include { MERGE_EXTENDED } from "../modules/merge_extended.nf"
include { BENCHMARKING } from "./benchmarking.nf"


integration_types = [
    "bbknn": "knn",
    "combat": "full",
    "desc": "embed",
    "harmony": "embed",
    "mnn": "full",
    "scanorama": "embed",
    "scanvi": "embed",
    "scvi": "embed",
    "trvaep": "embed",
    "scgen": "full",
    "unintegrated": "embed"
]

gpu_integrations = ["scgen"]

/**
 * Integrate individual datasets into a single-cell atlas
 *   - Perfom QC on individual datasets
 *   - Perform manual "seed" annotation of two datasets (one SS2, one 10x)
 *   - Perform data integration using SCANVI
 *   - Call doublets using solo
 */
workflow INTEGRATION {
    take:
        ch_adata
        ch_hvgs
        ch_integration_methods
        benchmark_hvgs

    main:
        ch_integration_methods = ch_integration_methods
            .filter{ it != "scvi" && it != "scanvi" }
            .branch{
                gpu: it in gpu_integrations
                cpu: it !in gpu_integrations
            }

        INTEGRATE(
            ch_adata,
            ch_hvgs,
            ch_integration_methods.cpu
        )

        INTEGRATE_GPU(
            ch_adata,
            ch_hvgs,
            ch_integration_methods.gpu
        )

        INTEGRATE_SCVI(
            ch_adata,
            ch_hvgs,
            "scvi"
        )

        ch_integrated = INTEGRATE.out.integrated
            .mix(INTEGRATE_GPU.out.integrated)
            .mix(INTEGRATE_SCVI.out.integrated)

        if (params.has_celltypes) {
            INTEGRATE_SCANVI(
                ch_adata,
                ch_hvgs,
                INTEGRATE_SCVI.out.model
            )

            ch_integrated = ch_integrated.mix(INTEGRATE_SCANVI.out.integrated)

            ch_model = INTEGRATE_SCANVI.out.model
            ch_scanvi_labels = INTEGRATE_SCANVI.out.labels
            ch_core_integrated = INTEGRATE_SCANVI.out.integrated
        } else {
            ch_model = INTEGRATE_SCVI.out.model
            ch_scanvi_labels = Channel.empty()
            ch_core_integrated = INTEGRATE_SCVI.out.integrated
        }

        ch_integrated_types = ch_integrated
            .map{ meta, adata -> [meta, adata, integration_types[meta.integration]] }

        BENCHMARKING(
            ch_adata,
            ch_integrated_types,
            benchmark_hvgs
        )

    emit:
        integrated = ch_integrated
        model = ch_model
        scanvi_labels = ch_scanvi_labels
}
