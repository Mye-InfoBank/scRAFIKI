include { INTEGRATE } from "../modules/integrate.nf"
include { INTEGRATE as INTEGRATE_GPU } from "../modules/integrate.nf"
include { INTEGRATE as INTEGRATE_SCVI } from "../modules/integrate.nf"
include { INTEGRATE_SCANVI } from "../modules/integrate_scanvi.nf"
include { EXTRACT_EMBEDDING } from "../modules/extract_embedding.nf"
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
        ch_preprocessed
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
            ch_preprocessed,
            ch_integration_methods.cpu
        )

        INTEGRATE_GPU(
            ch_preprocessed,
            ch_integration_methods.gpu
        )

        INTEGRATE_SCVI(
            ch_preprocessed,
            "scvi"
        )

        INTEGRATE_SCANVI(
            ch_preprocessed,
            INTEGRATE_SCVI.out.model
        )

        ch_integrated = INTEGRATE.out.integrated
            .mix(INTEGRATE_GPU.out.integrated)
            .mix(INTEGRATE_SCVI.out.integrated)
            .mix(INTEGRATE_SCANVI.out.integrated)

        ch_integrated_types = ch_integrated
            .map{ meta, adata -> [meta, adata, integration_types[meta.integration]] }

        EXTRACT_EMBEDDING(ch_integrated)

        BENCHMARKING(
            ch_preprocessed,
            ch_integrated_types,
            benchmark_hvgs
        )

    emit:
        integrated = ch_integrated
        scanvi_model = INTEGRATE_SCANVI.out.model
        obsm = EXTRACT_EMBEDDING.out
}
