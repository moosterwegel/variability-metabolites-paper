dir.create(here::here('data'))

path_covariates <- here::here('data', 'processed_covariate_data.csv')
path_lcms <- here::here('data', 'processed_lcms_data.csv')

path_annotations <- here::here('data', 'annotations.xlsx')
path_kegg <- here::here('data', 'annotations_plus_kegg_pathways.csv')
path_ancestors <- here::here('data', 'ancestors_annotations.xlsx')

curl::curl_download("https://zenodo.org/record/8156759/files/processed_covariate_data.csv?download=1", path_covariates)
curl::curl_download("https://zenodo.org/record/8156759/files/processed_lcms_data.csv?download=1", path_lcms)
curl::curl_download("https://zenodo.org/record/8156759/files/annotations.xlsx?download=1", path_annotations)
curl::curl_download("https://zenodo.org/record/8156759/files/annotations_plus_kegg_pathways.csv?download=1", path_kegg)
curl::curl_download("https://zenodo.org/record/8156759/files/ancestors_annotations.xlsx?download=1", path_ancestors)