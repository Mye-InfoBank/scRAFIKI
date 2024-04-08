workflow TRANSFER {
    take:
    ch_base
    ch_model
    ch_extension
    
    main:
    INTEGRATE_SCARCHES(
        ch_extension,
        ch_base.join(ch_model),
        params.has_celltypes,
        ch_hvgs
    )

    MERGE_EXTENDED(
        INTEGRATE_SCARCHES.out.integrated,
        ch_core_integrated
    )
    
    emit:
    
}