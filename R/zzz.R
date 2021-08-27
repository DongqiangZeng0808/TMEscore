



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://github.com/DongqiangZeng0808/", pkgname, "\n\n")

  citation <- paste0(" If you use ", pkgname, " in published research, please cite:\n",
                     " --------------------------------", "\n",
                     " Tumor microenvironment characterization in gastric cancer identifies prognostic and imunotherapeutically relevant gene signatures.","\n",
                     " Cancer Immunology Research, 2019, 7(5), 737-750", "\n",
                     " DOI: 10.1158/2326-6066.CIR-18-0436 ","\n" ,
                     " PMID: 30842092","\n" ,
                     " --------------------------------","\n",
                     " Tumor microenvironment evaluation promotes precise checkpoint immunotherapy of advanced gastric cancer.","\n",
                     " Journal for ImmunoTherapy of Cancer, 2021, 9(8), e002467","\n",
                     " DOI: 10.1136/jitc-2021-002467","\n",
                     " PMID: 34376552","\n",
                     " --------------------------------")

  packageStartupMessage(paste0(msg, citation))
}

