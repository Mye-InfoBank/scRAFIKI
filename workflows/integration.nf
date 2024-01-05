include { INTEGRATE } from "../modules/local/integrate.nf"
include { INTEGRATE as INTEGRATE_SCVI } from "../modules/local/integrate.nf"
include { INTEGRATE_SCANVI } from "../modules/local/integrate_scanvi.nf"


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
]

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


    main:
        ch_integration_methods = ch_integration_methods
            .filter{ it != "scvi" && it != "scanvi" }

        INTEGRATE(
            ch_preprocessed,
            ch_integration_methods
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
            .mix(INTEGRATE_SCVI.out.integrated)
            .mix(INTEGRATE_SCANVI.out.integrated)

        ch_integrated_types = ch_integrated
            .map{ meta, adata -> [meta, adata, integration_types[meta.integration]] }


    emit:
        integrated = ch_integrated
        integrated_types = ch_integrated_types
        scanvi_model = INTEGRATE_SCANVI.out.model
}
