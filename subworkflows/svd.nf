include { RIJK_ZSCORE   }   from    '../modules/local/rijk_zscore'
include { SVD_ZSCORE    }   from    '../modules/local/svd_zscore'

workflow SVD {
    take:
    dataname
    ch_pq 
    pin_S 
    pin_z 
    light 
    unfilt
    v2
    bound 
    temp
    svd_type
    
    main:
    // Step 1: Calculate RIJK z-score
    RIJK_ZSCORE (
        dataname,
        ch_pq,
        pin_S,
        pin_z,
        light,
        unfilt,
        v2,
        bound
    )

    // Step 2: Calculate SVD z-score
    SVD_ZSCORE (
        RIJK_ZSCORE.out.svd,
        temp,
        svd_type      
    )
    
    emit: 
    geneMats    = SVD_ZSCORE.out.geneMats
    svd         = SVD_ZSCORE.out.svd

}