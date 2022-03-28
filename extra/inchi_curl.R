inchi_curl <- 
  function(
           key,
           .id,
           type = "inchikey",
           get = "InChIkey",
           save = paste0(.id, ".csv")
           ){
    http = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/", type, "/")
    http_end = paste0("/property/", paste(get, collapse = ","), "/CSV > ")
    curl <- "curl -s --connect-timeout 20 --retry 100 --retry-delay 30 "
    curl_http <- paste0(curl, http)
    system(paste0(curl_http, key, http_end, save))
  }
int_inchi_curl <- 
  function(
           seq,
           type = "inchikey",
           get = "InChIkey",
           ...
           ){
    init <- dload[seq, "init"]
    .id <- dload[seq, ".id"]
    save <- paste0(.id, ".csv")
    while(class(try(read.csv(save), silent = T))[1] == "try-error"){
      inchi_curl(init, .id, 
                 get = get,
                 type = type,
                 ...)
    }
  }
