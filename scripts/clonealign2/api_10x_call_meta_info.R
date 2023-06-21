library(httr)

query_gsc_to_colossus_tenx <- function(str) {
  query <- paste0("colossus.molonc.ca/api/tenxlibrary/?name=", str)
  g <- GET(query, authenticate("htran", "3cbJNZHz"))
  json <- content(g)$results
  return(json)
}


meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata.csv')

View(meta)
sids <- c()
series <- c()
for(lid in meta$library_id){
  output <- query_gsc_to_colossus_tenx(lid)
  mouse_id <- output[[1]]$sample$sample_id
  if(is.null(mouse_id)){
    mouse_id <- ''
  }
  sids <- c(sids, mouse_id)
  series_id <- output[[1]]$sample$xenograft_id
  if(is.null(series_id)){
    series_id <- ''
  }
  series <- c(series, series_id)
}
df <- tibble::tibble(library_id=meta$library_id, sid_query=sids, series=series)
df
dim(df)
dim(meta)
df <- df[grepl('SA',df$sid_query),]
df
meta <- meta %>% left_join(df, by='library_id')
meta
sum(meta$mouse_id==meta$sid_query)
dim(meta)
meta$library_id[meta$mouse_id!=meta$sid_query]
meta$series
meta1 <- meta %>%
  dplyr::filter(mouse_id!=sid_query)
View(meta1)


data.table::fwrite(meta, '/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')




for(lid in meta1$library_id){
  output <- query_gsc_to_colossus_tenx(lid)
  
  
}

