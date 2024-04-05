workflow EXTEND {
    ch_adata = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obs = Channel.empty()
    
    emit:
        adata = ch_adata
        obsm = ch_obsm
        obs = ch_obs
}