#' Collapse genotype columns
#'
#' @param geno Integer genotype matrix with individuals in columns.
#' @return A character vector containing one genotype string per column.
#' @export
collapse_genotypes <- function(geno) {
  collapse_genotypes_cpp(geno)
}

#' Collapse genotype columns
#'
#' @param int1 int1
#' @param int2 int2
#' @return test
#' @export
kinship <- function(int1, int2) {
  kinship_cpp(int1, int2)
}

#' Find intervals application
#'
#' @param porc porc
#' @param parent_rec1 parent_rec1
#' @param parent_rec2 parent_rec2
#' @return test
#' @export

find_intervals_joint <- function(porc, parent_rec1, parent_rec2){
  find_intervals_joint_cpp(
    porc,
    parent_rec1,
    parent_rec2
  )
}

#' Breeding loops
#'
#' @param porc porc
#' @param start1 start1
#' @param start2 start2
#' @param end1 end1
#' @param end2 end2
#' @param parent_rec1 parent_rec1
#' @param parent_rec2 parent_rec2
#' @param parent_ori1 parent_ori1
#' @param parent_ori2 parent_ori2
#' @param mut_bool mut_bool
#' @param store_mut store_mut
#' @param parent_mut1 parent_mut1
#' @param parent_mut2 parent_mut2
#' @param recombination_precision recombination_precision
#' @param tmp1 tmp1
#' @return test
#' @export


breeding_loop <- function(porc,
                          start1, start2, end1, end2,
                          parent_rec1, parent_rec2, parent_ori1, parent_ori2,
                          mut_bool, store_mut,
                          parent_mut1, parent_mut2,
                          recombination_precision,
                          tmp1){

  breeding_loop_cpp(
    porcSEXP = porc,
    start1 = start1,
    start2 = start2,
    end1 = end1,
    end2 = end2,
    parent_rec1SEXP = parent_rec1,
    parent_rec2SEXP = parent_rec2,
    parent_ori1 = parent_ori1,
    parent_ori2 = parent_ori2,
    mut_bool = mut_bool,
    store_mut = store_mut,
    parent_mut1 = parent_mut1,
    parent_mut2 = parent_mut2,
    recombination_precision = recombination_precision,
    tmp1 = tmp1
  )
}

#' Breeding loops
#'
#' @param v1 v1
#' @param v2 v2
#' @param n1 n1
#' @return test
#' @export

merge_porc <- function(v1,v2,n1){
  merge_porc_cpp(
    v1,
    v2,
    n1
  )
}

