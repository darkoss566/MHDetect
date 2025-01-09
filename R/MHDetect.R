#' MHDetect Function
#'
#' This function detects and classifies indels (deletions, insertions, and MNVs)
#' as microhomology-mediated end joining (MMEJ)-dependent.
#'
#' @param vcf A VCF file (Variant Call Format) object read with readVcf.
#' @param k Number of base pairs to analyze before and after a deletion.
#' @param N Minimum number of matching nucleotides to classify as MMEJ.
#' @param genome The genome reference, e.g., BSgenome.Hsapiens.UCSC.hg19.
#' @param Interval Length of the region for analyzing DNA repair.
#' @return A data frame with classified indels and additional features.
#' @import GenomicRanges
#' @import Biostrings
#' @import BSgenome
#' @import VariantAnnotation
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' strelka <- "path/to/your/file.vcf.gz"
#' vcf_data <- readVcf(strelka, genome = "hg19")
#' result <- MHDetect(vcf_data, k = 25, N = 2,
#'                    genome = BSgenome.Hsapiens.UCSC.hg19, Interval = 25)
#' }
#' @export

MHDetect=function(vcf,k,N,genome, Interval){

  class_indel = function(vcf){
    # Get indel rows from VCF
    indels = rowRanges(vcf)

    # Filter indels with FILTER = 'PASS'
    indels_pass = indels[indels$FILTER == 'PASS']

    # Set sequence levels style to "UCSC"
    seqlevelsStyle(indels_pass) <- "UCSC"

    # Create a data frame from indels_pass
    indels_pass_data_frame <- data.frame(indels_pass)

    # Create data frames for indel properties
    Start <- data.frame(indels_pass_data_frame$start)
    colnames(Start) <- c('start')

    End <- data.frame(indels_pass_data_frame$end)
    colnames(End) <- c('end')

    ALT_indels = data.frame(indels_pass$ALT)[-c(1, 2)]
    REF_indels = data.frame(indels_pass$REF)

    chr <- data.frame(indels_pass_data_frame$seqnames)
    colnames(chr) <- c('seqnames')

    indels_data.frame = data.frame(chr, Start, End, ALT_indels, REF_indels)
    colnames(indels_data.frame) <- c('seqnames', 'start', 'end', 'ALT', 'REF')

    # Create a data frame for indel classification
    del_ins_MNV_SNV <- data.frame(matrix(nrow = nrow(indels_data.frame), ncol = 1))
    colnames(del_ins_MNV_SNV) <- c('indel')

    # Classify indels based on length and differences between ALT and REF
    for (i in 1:nrow(indels_data.frame)){
      length_ALT = nchar(indels_data.frame$ALT)[i]
      length_REF = nchar(indels_data.frame$REF)[i]

      if (length_REF > length_ALT){
        del_ins_MNV_SNV[i, ] = 'deletion'
      } else if (length_REF < length_ALT){
        del_ins_MNV_SNV[i, ] = 'insertion'
      } else if (indels_data.frame$ALT[i] != indels_data.frame$REF[i] && length_REF == 1){
        del_ins_MNV_SNV[i, ] = 'SNV'
      } else if (length_REF == length_ALT && indels_data.frame$ALT[i] == indels_data.frame$REF[i]){
        del_ins_MNV_SNV[i, ] = 'no_class'
      } else {
        del_ins_MNV_SNV[i, ] = 'MNV'
      }
    }

    indels_data.frame = data.frame(indels_data.frame, del_ins_MNV_SNV)

    # Return the classified indel data frame
    return(indels_data.frame)
  }

  #deletions
  deletions_maker = function(del_ALT){
    # Split ALT and REF for easier manipulation
    ALT_del = data.frame(del_ALT$ALT)
    colnames(ALT_del) <- c('ALT')

    REF_del = data.frame(del_ALT$REF)
    colnames(REF_del) <- c('REF')

    # Combine ALT and REF to create deletions
    deletions = data.frame(REF_del, ALT_del)

    # Create a data frame to store deletions
    deletions1 = data.frame(matrix(ncol = 1, nrow = nrow(deletions)))
    colnames(deletions1) <- c('deletion')

    # Loop to generate deletion sequences
    for (i in 1:nrow(deletions)){
      deletion = sub(deletions$ALT[i], '', deletions$REF[i])
      deletions1$deletion[i] = deletion
    }

    # Combine deletions1 with the original data frame
    deletions11 = data.frame(del_ALT, deletions1)

    colnames(deletions11) <- c("seqnames", "start", "end", "ALT", "REF", "indel_class", "deletion")

    return(deletions11)
  }
  # This code defines a function called 'deletions_before_after' that extracts sequences before and after deletions in genomic data.
  deletions_before_after = function(delecje_GRanges, k){
    # Keep only standard chromosomes
    Strelka_standard_chromosomes <- keepStandardChromosomes(delecje_GRanges, pruning.mode = 'coarse')
    delecje_dataframe <- data.frame(Strelka_standard_chromosomes)

    # Set the sequence style
    seqlevelsStyle(Hsapiens) <- "UCSC"

    # Deletion
    now = GRanges(delecje_dataframe$seqnames, IRanges(start = delecje_dataframe$start + 1, end = delecje_dataframe$end))
    seqlevelsStyle(now) <- "UCSC"
    seq_now = data.frame(getSeq(Hsapiens, now))
    colnames(seq_now) <- c('deletion')

    # Sequence after the deletion
    after = GRanges(delecje_dataframe$seqnames, IRanges(start = delecje_dataframe$end + 1, end = delecje_dataframe$end + k))
    seqlevelsStyle(after) <- "UCSC"
    seq_after = data.frame(getSeq(Hsapiens, after))
    colnames(seq_after) <- c('after')

    # Sequence before the deletion
    before = GRanges(delecje_dataframe$seqnames, IRanges(start = delecje_dataframe$start + 1 - k, end = delecje_dataframe$start))
    seqlevelsStyle(before) <- "UCSC"
    seq_before = data.frame(getSeq(Hsapiens, before))
    colnames(seq_before) <- c('before')

    # Create a new data frame
    new_data_frame = data.frame(Strelka_standard_chromosomes, seq_before, seq_now, seq_after)

    return(new_data_frame)
  }

  #insertions
  insertions_maker = function(ins_ALT) {
    # Divide ALT and REF for easier manipulation
    ALT_ins = data.frame(ins_ALT$ALT)
    colnames(ALT_ins) <- c('ALT')
    REF_ins = data.frame(ins_ALT$REF)
    colnames(REF_ins) <- c('REF')
    insertions = data.frame(REF_ins, ALT_ins)

    # Create a data frame to store insertions
    insertions1 = data.frame(matrix(ncol = 1, nrow = nrow(insertions)))
    colnames(insertions1) <- c('insercja')

    # Loop to create insertion sequences
    for (i in 1:nrow(insertions)) {
      insertion = sub(insertions$REF[i], '', insertions$ALT[i])
      insertions1$insercja[i] = insertion
    }

    insertions11 = data.frame(ins_ALT, insertions1$insercja)

    colnames(insertions11) <- c("seqnames", "start", "end", "ALT", "REF", "indel_class", "indels")

    return(insertions11)
  }
  # This code defines a function called 'insertions_before_after' that extracts sequences before and after insertions in genomic data.
  insertions_before_after=function(insercje_GRanges, k){

    insercje_GRanges<-GRanges_insercje

    # Filtering and keeping only standard chromosomes in the genomic data
    Strelka_standard_chromosomes <- keepStandardChromosomes(insercje_GRanges, pruning.mode= 'coarse')

    # Creating a data frame from the filtered genomic data
    insercje_dataframe <- data.frame(Strelka_standard_chromosomes)

    # Setting the style of sequence levels for the human genome to "UCSC"
    seqlevelsStyle(Hsapiens)<-"UCSC"

    # Generating sequences after each insertion position
    now = GRanges(insercje_dataframe$seqnames, IRanges(start=insercje_dataframe$start+1, end=insercje_dataframe$end))
    seqlevelsStyle(now)<-"UCSC"

    # Extracting sequences using the human genome reference
    seq_now = data.frame(getSeq(Hsapiens, now))
    colnames(seq_now) <- c('deletion')

    # Generating sequences after each insertion position (specified by 'k' length)
    after = GRanges(insercje_dataframe$seqnames, IRanges(start=insercje_dataframe$end+1, end=insercje_dataframe$end+k))
    seqlevelsStyle(after)<-"UCSC"

    # Extracting sequences using the human genome reference
    seq_after = data.frame(getSeq(Hsapiens, after))
    colnames(seq_after) <- c('after')

    # Generating sequences before each insertion position (specified by 'k' length)
    before = GRanges(insercje_dataframe$seqnames, IRanges(start=insercje_dataframe$start+1-k, end=insercje_dataframe$start))
    seqlevelsStyle(before)<-"UCSC"

    # Extracting sequences using the human genome reference
    seq_before = data.frame(getSeq(Hsapiens, before))
    colnames(seq_before) <- c('before')

    # Creating a final data frame that includes relevant information
    new_data_frame=data.frame(Strelka_standard_chromosomes, seq_before, insercje$indels, seq_after)

    # Returning the final data frame
    return(new_data_frame)

  }

  #MNV
  #Filters MNV data, extracts sequences before and after each position, and compiles relevant information into a data frame.
  mnv_before_after=function(mnv_GRanges, k){

    #mnv_GRanges<-GRanges_mnv

    # Filtering and keeping only standard chromosomes in the genomic data
    Strelka_standard_chromosomes <- keepStandardChromosomes(mnv_GRanges, pruning.mode= 'coarse')

    # Creating a data frame from the filtered genomic data
    mnv_dataframe <- data.frame(Strelka_standard_chromosomes)

    # Setting the style of sequence levels for the human genome to "UCSC"
    seqlevelsStyle(Hsapiens)<-"UCSC"

    # Generating sequences after each mnv position
    now = GRanges(mnv_dataframe$seqnames, IRanges(start=mnv_dataframe$start+1, end=mnv_dataframe$end))
    seqlevelsStyle(now)<-"UCSC"

    # Extracting sequences using the human genome reference
    seq_now = data.frame(getSeq(Hsapiens, now))
    colnames(seq_now) <- c('deletion')

    # Generating sequences after each mnv position (specified by 'k' length)
    after = GRanges(mnv_dataframe$seqnames, IRanges(start=mnv_dataframe$end+1, end=mnv_dataframe$end+k))
    seqlevelsStyle(after)<-"UCSC"

    # Extracting sequences using the human genome reference
    seq_after = data.frame(getSeq(Hsapiens, after))
    colnames(seq_after) <- c('after')

    # Generating sequences before each insertion position (specified by 'k' length)
    before = GRanges(mnv_dataframe$seqnames, IRanges(start=mnv_dataframe$start+1-k, end=mnv_dataframe$start))
    seqlevelsStyle(before)<-"UCSC"

    # Extracting sequences using the human genome reference
    seq_before = data.frame(getSeq(Hsapiens, before))
    colnames(seq_before) <- c('before')

    # Creating a final data frame that includes relevant information
    new_data_frame=data.frame(Strelka_standard_chromosomes, seq_before, mnv_ALT$ALT, seq_after)


    # Returning the final data frame
    return(new_data_frame)

  }
  #Generates a data with position of DSB around MNV +/- k, between each nucleotides.
  create_DSB<-function(specified_range){

    specified_range2<-data.frame(specified_range)
    dsb_data_frame<-data.frame()

    for (j in 1:nrow(specified_range2)){
      mnv_row<-specified_range2[j,]
      seq<-mnv_row$seqnames
      indel<-mnv_row$indels

      for (i in (mnv_row$start-1):(mnv_row$end+1)){
        dsb<-GRanges(seqnames=seq, IRanges(start=i, end=i))
        dsb_data_frame<-rbind(dsb_data_frame, data.frame(dsb, indel))
      }
    }

    return(dsb_data_frame)

  }
  rangeMid = function(GRangesIn,N=10) {
    mid = round((start(GRangesIn)+end(GRangesIn))/2)
    outGR = GRangesIn
    start(outGR) = mid - N/2
    end(outGR) = mid + N/2
    return(outGR)
  }
  coordinates_overlap = function(s1,e1,s2,e2) {
    return(ifelse((s1<s2 & e1>=s2) | (s1>s2 & e2>=s1) | s1==s2, TRUE, FALSE))
  }
  count_common_nucleotides = function(seq1,seq2) {
    if (nchar(seq1)==0 | nchar(seq2)==0) {
      ct = NA
    } else {
      ct = 0
      for (i in 1:min(nchar(seq1),nchar(seq2))) {
        if (substr(seq1,i,i)==substr(seq2,i,i)) {
          ct=ct+1
        } else {
          break;
        }
      }
    }
    return(ct)
  }
  mhAlignment = function(randSeq5fill,randSeq3fill) {

    randSeq5fill = Biostrings::DNAString(randSeq5fill)
    randSeq3fill = Biostrings::DNAString(randSeq3fill)

    #check the length of the sequences
    SeqLengths = c(randSeq5fill@length,randSeq3fill@length)
    MaxOffset = min(SeqLengths) -1
    LenDiff = SeqLengths[1]-SeqLengths[2]


    #offset starting position due to uneven length of the sequences
    if (SeqLengths[1]>SeqLengths[2]) {
      start_offset = LenDiff
    } else {
      start_offset = 0
    }

    #identify matching bases
    s5_pos_end = s5_pos_start = s3_pos_end = s3_pos_start = best_mh = NA
    best_match = 0
    #for each offset of the sequence (move the randSeq3fill by 1bp in relation to randSeq5fill)
    #  the loop iterates the elements in reverse to prefer binding in the vicinity of DSB
    for (ofs in MaxOffset:1) {
      s_match = 0
      s_matchpos = NA
      #check each base of the two sequences in places where there is a sequence to compare
      for (step in 1:(MaxOffset-ofs+1)) {
        if (randSeq5fill[start_offset+ofs+step] == randSeq3fill[step]) {
          s_match = s_match + 1
          #save the starting position if we have a match
          if (s_match==1) {
            s_matchpos = step
          }
          #save the best result
          if (s_match>best_match) {
            best_match = s_match
            s5_pos_end = start_offset+ofs+step

            s5_pos_start = s5_pos_end-s_match+1
            s3_pos_end = step
            s3_pos_start = s3_pos_end-s_match+1
            best_mh = randSeq5fill[s5_pos_start:s5_pos_end]
          }
        } else {
          s_match = 0
        }
      }
    }
    return(list(s5_pos_start = s5_pos_start, s5_pos_end = s5_pos_end, s3_pos_start = s3_pos_start, s3_pos_end = s3_pos_end, mh_seq = as.character(best_mh)))
  }
  #Simulates the template-switching mechanism in DNA repair.
  simulate_template_switch = function(DSB_create) {

    #generate random positions with DSB
    randPos<-makeGRangesFromDataFrame(DSB_create, keep.extra.columns=TRUE)
    randPosMid = rangeMid(randPos,Interval*2)
    randSeq = Biostrings::getSeq(Hsapiens,randPosMid)

    #5’-3’ end-resection
    randPos$randSeq5 = subseq(randSeq,1,Interval)
    randPos$randSeq3 = subseq(randSeq,Interval+1,Interval*2)


    #randPos$repair_strand = sample(c('F','R'),length(randPos),replace = F)

    col_names<-colnames(data.frame(randPos))


    randPos_F<-data.frame(randPos,'F')
    colnames(randPos_F)<-append(col_names,'repair_strand')

    randPos_R<-data.frame(randPos,'R')
    colnames(randPos_R)<-append(col_names,'repair_strand')

    data<-rbind(randPos_F, randPos_R)
    randPos=makeGRangesFromDataFrame(data,keep.extra.columns=TRUE)

    l<-list('start_seq','end_seq')
    col_names<-append(l,col_names)

    randPos_F<-data.frame(start(randPosMid),end(randPosMid),randPos,'F')
    colnames(randPos_F)<-append(col_names,'repair_strand')

    randPos_R<-data.frame(start(randPosMid),end(randPosMid),randPos,'R')
    colnames(randPos_R)<-append(col_names,'repair_strand')

    data2<-rbind(randPos_F, randPos_R)


    ReportTable = data.frame()
    for (i in 1:length(randPos)) {

      #i=871
      randSeq5 = randPos$randSeq5[[i]]
      randSeq3 = randPos$randSeq3[[i]]
      repair_strand = randPos$repair_strand[i]

      OriginalSequence = c(randSeq5,randSeq3)
      OriginalSequence <- paste(OriginalSequence, collapse = "")
      OriginalSequence = DNAString(OriginalSequence)

      ## 1. Anneling of microhomologies
      mh1_pos = mhAlignment(randSeq5,randSeq3)  #F/R selection doesn't matter here since both sequences have equal length

      mh1seq = mh1_pos$mh_seq
      mh1len = nchar(mh1seq)

      if (!is.na(mh1_pos$mh_seq)){

        ## 2. Flap removal
        randSeq5rm = subseq(randSeq5,1,mh1_pos$s5_pos_end)
        randSeq3rm = subseq(randSeq3,mh1_pos$s3_pos_start,Interval)
        #removed sequence fragments
        randSeq5rm_RM = subseq(randSeq5,mh1_pos$s5_pos_end+1,Interval)
        randSeq3rm_RM = subseq(randSeq3,1,mh1_pos$s3_pos_start-1)

        ## 3. Fill-in synthesis (only one, random direction)
        #pick strand
        filled_seq = Biostrings::DNAString(paste0(randSeq5rm,substring(Biostrings::DNAStringSet(randSeq3rm),mh1len+1)))
        if (repair_strand=="F") {
          randSeq5fill = filled_seq  #this is the same for F and R since we do not use complement sequence
          randSeq3fill = randSeq3rm
        } else {
          randSeq5fill = randSeq5rm
          randSeq3fill = filled_seq
        }

        ## 4. new MH (here pairwiseAlignment wont work since the best will always use the previous MH)
        mh2_pos = mhAlignment(randSeq5fill,randSeq3fill)

        mh2seq = mh2_pos$mh_seq
        mh2len = nchar(mh2seq)

        if (!is.na(mh2_pos$mh_seq)){

          ## 5. 2nd flap removal
          randSeq5rm2 = subseq(randSeq5fill,1,mh2_pos$s5_pos_end)
          randSeq3rm2 = subseq(randSeq3fill,mh2_pos$s3_pos_start,nchar(randSeq3fill))

          ## 6. 2nd fill-in synthesis  (both directions) - the sequences are now identical
          randSeq5fill2 = Biostrings::DNAString(paste0(randSeq5rm2,substring(as.character(randSeq3rm2),mh2len+1)))
          randSeq3fill2 = randSeq5fill2 # or DNAString(paste0(substring(as.character(randSeq5rm2),1,mh2_pos$s5_pos_start-1),randSeq3rm2))
          FinalSequence = randSeq5fill2

          ## 7. determine the impact of both MH repairs
          mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -10, baseOnly = TRUE)
          aln = pairwiseAlignment(OriginalSequence,FinalSequence,type='global',substitutionMatrix = mat,gapOpening=10,gapExtension=1)


          ## statistics:
          #length of the flap removal sequences from the first step + MH1 length -OR- the distance between MH1 motifs in the sequence NNNNNN[MH1]NNNNNNNN[MH1]NNNNNN + MH1 length
          deletion_seq = c(Biostrings::DNAString(mh1seq),Biostrings::DNAString(randSeq5rm_RM),Biostrings::DNAString(randSeq3rm_RM))
          deletion_len = length(deletion_seq)
          if(mh1_pos$s5_pos_start+mh2_pos$s3_pos_start-1<=mh2_pos$s5_pos_start-1){
            #MH2 len + fill in from F or R  -OR-  the distance between MH1 and MH2 motifs in the sequence NNNNNN[MH1]NNNNNNNN[MH2]NNNNNN
            if (repair_strand=='F') {
              #insertion_seq = subseq(randSeq5fill2,mh1_pos$s5_pos_end+1,mh2_pos$s5_pos_start-1)
              insertion_seq = subseq(randSeq5fill2,mh1_pos$s5_pos_start+mh2_pos$s3_pos_start-1,mh2_pos$s5_pos_start-1)
            } else {
              #insertion_seq = subseq(randSeq3fill2,mh2_pos$s3_pos_end+1,mh2_pos$s5_pos_end)
              insertion_seq = subseq(randSeq3fill2,mh2_pos$s3_pos_start,mh2_pos$s5_pos_start-1)
            }

            insertion_len = length(insertion_seq)
            #original len
            original_len = length(OriginalSequence)
            #final len
            observed_len = length(randSeq5fill2)

          }else{

            insertion_seq=0
            insertion_len = 0
            #original len
            original_len = length(OriginalSequence)
            #final len
            observed_len = length(randSeq5fill2)

          }

          #final_len == original_len - deletion_len + insertion_len

          tReportTable = data.frame(DSB_chr = as.character(GenomeInfoDb::seqnames(randPos[i])), #chromosome in which DSB occurs
                                    DSB_pos = start(randPos[i]),                                 #position at which DSB occurs
                                    start_seq = data2$start_seq[i],
                                    end_seq = data2$end_seq[i],
                                    indel = data2$indel[i],
                                    original_context_5p = as.character(randSeq5),                #context of the DSB located upstream (towards 5' end)
                                    original_context_3p = as.character(randSeq3),                #context of the DSB located downstream (towards 3' end)
                                    original_seq = as.character(OriginalSequence),               #combined context sequences upstream + downstream
                                    observed_seq = as.character(FinalSequence),                  #sequence obtained as a result of repair involving template switching
                                    deletion_seq = as.character(deletion_seq),                   #sequence of the introduced deletion
                                    insertion_seq = as.character(insertion_seq),                 #sequence of the introduced insertion
                                    observed_len,                                                #length of the observed_seq
                                    deletion_len,                                                #length of the deletion
                                    insertion_len,                                               #length of the insertion
                                    first_fill_in_strand = repair_strand,                        #strand at which the repair first occurred before template switch (F/R)
                                    mh1seq,                                                      #sequence of the first microhomology motif (MH1)
                                    mh2seq,                                                      #sequence of the second microhomology motif (MH2)
                                    mh1_start1 = mh1_pos$s5_pos_start,                           #start position of the first occurrence of MH1 inside the original_seq
                                    mh1_end1 = mh1_pos$s5_pos_end,                               #end position of the first occurrence of MH1 inside the original_seq
                                    mh1_start2 = Interval + mh1_pos$s5_pos_start,                #start position of the second occurrence of MH1 inside the original_seq
                                    mh1_end2 = Interval + mh1_pos$s5_pos_end,                    #end position of the second occurrence of MH1 inside the original_seq

                                    mh2_start1 = ifelse(repair_strand=='F',                      #start position of the first occurrence of MH2 inside the original_seq
                                                        deletion_len + mh2_pos$s5_pos_start,
                                                        mh2_pos$s5_pos_start),
                                    mh2_end1 = ifelse(repair_strand=='F',                        #end position of the first occurrence of MH2 inside the original_seq
                                                      deletion_len + mh2_pos$s5_pos_end,
                                                      mh2_pos$s5_pos_end),
                                    mh2_start2 = ifelse(repair_strand=='F',                      #start position of the second occurrence of MH2 inside the original_seq
                                                        length(randSeq3rm_RM) + Interval + mh2_pos$s3_pos_start,
                                                        mh2_pos$s3_pos_start),
                                    mh2_end2 = ifelse(repair_strand=='F',                        #end position of the second occurrence of MH2 inside the original_seq
                                                      length(randSeq3rm_RM) + Interval + mh2_pos$s3_pos_end,
                                                      mh2_pos$s3_pos_end),

                                    final_alignment_original = as.character(aln@pattern),        #result of global alignment between original_seq and observed_seq: original_seq part
                                    final_alignment_observed = as.character(aln@subject))        #result of global alignment between original_seq and observed_seq: observed_seq part

          ReportTable = rbind(ReportTable,tReportTable)


        }

      }

    }




    #additional statistics for the simulated sites
    ReportTable_ext = ReportTable %>% rowwise() %>%
      #MH overlaps
      mutate(mh1_overlap = coordinates_overlap(mh1_start1,mh1_end1,mh1_start2,mh1_end2),             #True if the MH1 pairs overlap (this should never occur)
             mh2_overlap = coordinates_overlap(mh2_start1,mh2_end1,mh2_start2,mh2_end2),             #True if the MH2 pairs overlap
             mh1_mh2_overlap_P = coordinates_overlap(mh1_start1,mh1_end1,mh2_start1,mh2_end1),       #True if the MH1 and MH2 pairs overlap - first motivs from pairs
             mh1_mh2_overlap_M = coordinates_overlap(mh1_start2,mh1_end2,mh2_start2,mh2_end2),       #True if the MH1 and MH2 pairs overlap - second motivs from pairs
             mh1_mh2_overlap_PM = coordinates_overlap(mh1_start2,mh1_end2,mh2_start2,mh2_end2)) %>%  #True if the MH1 and MH2 pairs overlap - first and second motiv from pairs
      mutate(mh1_mh2_overlap = mh1_mh2_overlap_P | mh1_mh2_overlap_M | mh1_mh2_overlap_PM) %>%  #True if any of the  MH1 and MH2 motives overlap

      #alignment details
      mutate(observed_insLen = str_count(final_alignment_original, "-"),        #number of gaps in the original sequence (sum of all inserted nucleotides)
             observed_delLen = str_count(final_alignment_observed, "-")) %>%    #number of gaps in the observed sequence (sum of all deleted nucleotides)
      mutate(observed_exp_ins_match = observed_insLen==insertion_len,     #True if the sum of all inserted nucleotides (from alignment) is equal to observed len of insertion
             observed_exp_del_match = observed_delLen==deletion_len) %>%  #True if the sum of all deleted nucleotides (from alignment) is equal to observed len of deletion

      #position of insertion and deletion
      mutate(deletion_start = DSB_pos + mh1_start1 - Interval - 1,                   #start position of the deletion in reference genome
             deletion_end = DSB_pos + mh1_start1 + deletion_len - Interval - 2) %>%  #end position of the deletion in reference genome
      mutate(insertion_start = ifelse(first_fill_in_strand=='F',                     #start position of the insertion in reference genome
                                      DSB_pos + mh2_start2 - Interval - 2,
                                      DSB_pos + mh2_start1 - Interval - 2),
             insertion_end = ifelse(first_fill_in_strand=='F',                       #end position of the insertion in reference genome
                                    DSB_pos + mh2_start2 - Interval- 1,
                                    DSB_pos + mh2_start1 - Interval- 1)) %>%

      ## sequence validation
      mutate(deletion_seq_control=as.character(Biostrings::getSeq(Hsapiens,GenomicRanges::GRanges(DSB_chr,IRanges::IRanges(deletion_start,deletion_end)))[[1]])) %>%  #extracted deletion sequence based on calculated coordinates
      mutate(deletion_seq_match = deletion_seq_control==deletion_seq) %>%
      mutate(insertion_seq_context=paste(as.character(Biostrings::getSeq(Hsapiens,GenomicRanges::GRanges(DSB_chr,IRanges::IRanges(insertion_start-2,insertion_end-1)))[[1]]),
                                         as.character(Biostrings::getSeq(Hsapiens,GenomicRanges::GRanges(DSB_chr,IRanges::IRanges(insertion_start+1,insertion_end+2)))[[1]]),sep='_'))  %>%

      ## detection
      # del_up_mh sequence of mh
      mutate(del_up_mh = as.character(Biostrings::getSeq(Hsapiens,GenomicRanges::GRanges(DSB_chr,IRanges::IRanges(deletion_end+1,deletion_end+nchar(mh1seq))))[[1]]),
             ins_up_mh = as.character(Biostrings::getSeq(Hsapiens,GenomicRanges::GRanges(DSB_chr,IRanges::IRanges(insertion_end,insertion_end+nchar(mh2seq)-1)))[[1]])) %>%
      # sequence agreement
      mutate(mh1_is_up_del = del_up_mh == mh1seq,
             mh2_is_up_ins = ins_up_mh == mh2seq,
             mh1_match_up_del = count_common_nucleotides(del_up_mh,mh1seq),
             mh2_match_up_ins = count_common_nucleotides(ins_up_mh,mh2seq),
             del_match_up_del = count_common_nucleotides(del_up_mh,deletion_seq),
             ins_match_up_ins = count_common_nucleotides(ins_up_mh,insertion_seq)) %>%

      #pattern is detectable
      mutate(is_detectable = del_match_up_del>1 & ins_match_up_ins>1)
    return(ReportTable_ext)
  }
  #Compares MNV and template-switched MNV sequences, calculating the percentage of agreement.
  FilterFunction <- function(ReportTable_data, mnv_before_after, k) {

    # Initialize an empty object ReportTable_grouped
    ReportTable_grouped <- c()
    # Loop through all rows of the mnv_before_after dataframe
    for (i in 1:nrow(mnv_before_after)) {
      # Select the i-th row from the mnv_before_after dataframe
      mnv_i <- mnv_before_after[i, ]
      #mnv_i <- mnv_before_after[which(mnv_before_after$indels=='TTCATA'),]

      # Select rows from the ReportTable_data where indel equals indels from mnv_i
      ReportTable_data_indel <-
        ReportTable_data[which(ReportTable_data$indel == mnv_i$indels), ]

      # Filter data within the indel group
      ReportTable_data_grouped <- ReportTable_data_indel[which(ReportTable_data_indel$start_seq >= (mnv_i$start) - k & ReportTable_data_indel$end_seq <= (mnv_i$end) + k),]

      # Add results to ReportTable_grouped
      ReportTable_grouped <-
        rbind(ReportTable_grouped, ReportTable_data_grouped)
    }


    # Initialize an empty object ReportTable_big_data
    ReportTable_big_data <- c()

    # Loop through all rows of the ReportTable_grouped dataframe
    for (i in 1:nrow(ReportTable_grouped)) {

      # Select the i-th row from the ReportTable_grouped dataframe
      indel_chose <- ReportTable_grouped[i, ]

      # Select rows from the mnv_before_after dataframe satisfying the conditions
      #mnv_before_after$indels == indel_chose$indel: This condition ensures that
      #the indels column in mnv_before_after matches the indel value from the current row in indel_chose.
      #mnv_before_after$seqnames == indel_chose$DSB_chr: This condition ensures that
      #the seqnames column in mnv_before_after matches the DSB_chr value from the current row in indel_chose.
      #mnv_before_after$start >= (indel_chose$start_seq - k) & mnv_before_after$start <= (indel_chose$start_seq + k):
      #These conditions define a range for the start column in mnv_before_after.
      #Rows are included if their start value falls within the range defined by (indel_chose$start_seq - k) and
      #(indel_chose$start_seq + k)

      mnv_o <- mnv_before_after[which(
        mnv_before_after$indels == indel_chose$indel & #This condition ensures that the indels column in mnv_before_after matches the indel value from the current row in indel_chose.
          mnv_before_after$seqnames == indel_chose$DSB_chr &  #This condition ensures that the seqnames column in mnv_before_after matches the DSB_chr value from the current row in indel_chose.
          mnv_before_after$start >= (indel_chose$start_seq - k) &
          mnv_before_after$start <= (indel_chose$start_seq + k)), ] #These conditions define a range for the start column in mnv_before_after.
      #Rows are included if their start value falls within the range defined by (indel_chose$start_seq - k) and
      #(indel_chose$start_seq + k)

      if (nrow(mnv_o)>0){

        # Create indel_from_template_switching sequence
        indel_from_template_switching <-substr(
          indel_chose$observed_seq, #This is the observed sequence from the indel_chose dataframe. It represents the sequence where the indel event is observed.
          mnv_o$start - indel_chose$start_seq + 1, #This calculates the starting position within the observed sequence for the indel from the template switching. It adjusts the start position relative to the start_seq of the indel in indel_chose
          mnv_o$start - indel_chose$start_seq + nchar(indel_chose$indel)) #This calculates the ending position within the observed sequence for the indel from the template switching. It adjusts the end position relative to the start_seq of the indel in indel_chose and considers the length of the indel.
        str_indel_from_template_switching <-strsplit(indel_from_template_switching, "")[[1]]
        str_indel_chose <- strsplit(indel_chose$indel, "")[[1]]

        # Compare letters one by one and calculate the percentage of agreement
        suppressWarnings({
          vector <- str_indel_chose == str_indel_from_template_switching
          procent_TRUE <- sum(vector) / nchar(indel_chose$indel) * 100
        })

        # Create a dataframe and add it to ReportTable_big_data
        data_frame <-data.frame(indel_chose$DSB_pos,indel_chose$first_fill_in_strand,mnv_o, indel_from_template_switching, procent_TRUE)
        ReportTable_big_data <- rbind(ReportTable_big_data, data_frame)

      }


    }

    colnames(ReportTable_big_data)[1]<-'DSB_pos'
    colnames(ReportTable_big_data)[2]<-'first_fill_in_strand'
    List<-list(ReportTable_grouped,ReportTable_big_data)

    return(List)
  }




  # This function, 'REPEATS_counter', counts the number of repeats in sequences before and after insertions.
  REPEATS_counter=function(indels_before_after){

    # Create an empty data frame to store repeat counts
    powtorzenia = data.frame(matrix(ncol = 1, nrow = 0))
    colnames(powtorzenia) <- c('repeat_count')
    #microhomology = data.frame(matrix(ncol = 1, nrow = 0))
    #colnames(microhomology) <- c('microhomology')

    # Loop through each insertion sequence
    for (i in 1:nrow(indels_before_after)){

      length_indels<- nchar(toString(indels_before_after$indels[i]))

      if (length_indels > 0){

        # Limit the deletion length to the length of the 'after' sequence
        length_after <- nchar(indels_before_after$after[i])
        df <- substring(indels_before_after$indels[i], 1, length_after)

        # Find positions of repeats
        repeats = unlist(gregexpr(df, indels_before_after$after[i]))

        if (repeats[1] == 1){ # Check if the first position in the 'repeats' vector is 1, allowing for further analysis

          length_indels = nchar(toString(indels_before_after$indels[i]))
          a = diff(repeats) == length_indels # Calculate the difference between elements
          b = match(FALSE, a) # Count the number of repeats until FALSE

          if (is.na(b) == TRUE) {

            repeat_count = length(repeats) # Count of repeats from 'repeats'

            if (repeat_count >= 2){
              powtorzenia[nrow(powtorzenia) + 1, ] = repeat_count
              #microhomology[nrow(microhomology) + 1, ] = 'Repeat-mediated'
            } else{
              powtorzenia[nrow(powtorzenia) + 1, ] = repeat_count
              #microhomology[nrow(microhomology) + 1, ] = 'Another classification'
            }


          } else {

            repeat_count = b # Count of repeats until the first FALSE

            if (repeat_count >= 2){
              powtorzenia[nrow(powtorzenia) + 1, ] = repeat_count
              #microhomology[nrow(microhomology) + 1, ] = 'Repeat-mediated'
            } else{
              powtorzenia[nrow(powtorzenia) + 1, ] = repeat_count
              #microhomology[nrow(microhomology) + 1, ] = 'Another classification'
            }
          }


        } else {

          powtorzenia[nrow(powtorzenia) + 1, ] = 0
          #microhomology[nrow(microhomology) + 1, ] = 'Another classification'

        }

      } else {
        powtorzenia[nrow(powtorzenia) + 1, ] = 0
      }

    }

    new_data_frame = data.frame(indels_before_after, powtorzenia) #,microhomology)
    return(new_data_frame)
  }

  # This function, 'NUCLEOTIDE_counter', counts the number of matching nucleotides in sequences before and after insertions.
  NUCLEOTIDE_counter=function(REPEATS){

    # Create an empty data frame to store nucleotide match counts
    l_nukleotydow = data.frame(matrix(ncol = 1, nrow = 0))
    colnames(l_nukleotydow) <- c('matching_nucleotides_count')
    #n_klasyfikacja = data.frame(matrix(ncol = 1, nrow = 0))
    #colnames(n_klasyfikacja) <- c('n_classification')

    # Loop through each insertion sequence
    for (i in 1:nrow(REPEATS)){

      # Limit the length of deletion and after sequences to their respective lengths
      length_after <- nchar(REPEATS$after[i])
      length_indel <- nchar(REPEATS$indels[i])
      delecja_skrocona <- substring(REPEATS$indels[i], 1, length_after)
      after_skrocona <- substring(REPEATS$after[i], 1, length_indel)

      # Convert characters to data frames using 'utf8ToInt' and 'intToUtf8'
      indel <- as.data.frame(Map(intToUtf8, as.list(utf8ToInt(delecja_skrocona))))
      after <- as.data.frame(Map(intToUtf8, as.list(utf8ToInt(after_skrocona))))
      dl_del <- length(indel)

      # Compare each character in 'indel' and 'after' sequences
      line <- indel == after
      a <- match(FALSE, line)

      if (is.na(a) == TRUE){
        ile_true <- length(line)

        if (ile_true > REPEATS$repeat_count[i]){
          l_nukleotydow[nrow(l_nukleotydow) + 1, ] = ile_true
          #n_klasyfikacja[nrow(n_klasyfikacja) + 1,]  = 'Microhomology mediated'
        } else{
          l_nukleotydow[nrow(l_nukleotydow) + 1, ] = ile_true
          #n_klasyfikacja[nrow(n_klasyfikacja) + 1,]  = '0'
        }

      } else {
        ile_true <- a - 1

        if (ile_true > REPEATS$repeat_count[i]){
          l_nukleotydow[nrow(l_nukleotydow) + 1, ] = ile_true
          #n_klasyfikacja[nrow(n_klasyfikacja) + 1,]= 'Microhomology mediated'
        } else{
          l_nukleotydow[nrow(l_nukleotydow) + 1, ] = ile_true
          #n_klasyfikacja[nrow(n_klasyfikacja) + 1,]  = '0'

        }
      }

    }

    new_data_frame = data.frame(REPEATS, l_nukleotydow) #, n_klasyfikacja)
    return(new_data_frame)
  }

  # This function, 'classification', classifies insertion sequences based on the number of matching nucleotides and repeat counts.
  classification=function(nucleotide, N){

    # Create an empty data frame to store microhomology classification
    microhomology = data.frame(matrix(ncol = 1, nrow = (nrow(nucleotide))))
    colnames(microhomology) <- c('microhomology_classification')

    # Loop through each insertion sequence
    for (i in 1:nrow(nucleotide)){

      # Check if the repeat count is greater than or equal to 2
      if ((nucleotide$repeat_count[i]) >= 2){
        microhomology[i,] = 'Repeat-mediated'

      }else{

        # Check if the count of matching nucleotides is greater than or equal to N
        if ((nucleotide$matching_nucleotides_count[i]) >= N){
          microhomology[i,] = 'Microhomology mediated'

        }else{
          microhomology[i,] = 'No classification'

        }

      }
    }

    # Combine the original nucleotide data and microhomology classification
    new_data_frame = data.frame(nucleotide, microhomology)
    return(new_data_frame)
  }

  # Extract indel data from the VCF file
  indels_data.frame <- class_indel(vcf)

  # Separate deletions (DEL) and insertions (INS) based on 'indel' column
  del_ALT <- if (nrow(indels_data.frame[indels_data.frame$indel == 'deletion',]) > 0) {
    indels_data.frame[indels_data.frame$indel == 'deletion',]
  }

  ins_ALT <- if (nrow(indels_data.frame[indels_data.frame$indel == 'insertion',]) > 0) {
    indels_data.frame[indels_data.frame$indel == 'insertion',]
  }

  snv_ALT <- if (nrow(indels_data.frame[indels_data.frame$indel == 'SNV',]) > 0) {
    indels_data.frame[indels_data.frame$indel == 'SNV',]
  }

  mnv_ALT <- if (nrow(indels_data.frame[indels_data.frame$indel == 'MNV',]) > 0) {
    indels_data.frame[indels_data.frame$indel == 'MNV',]
  }


  # Process deletions (DEL) if there are any
  if (!is.null(del_ALT) && nrow(del_ALT) != 0){
    # Create deletion data
    delecje <- deletions_maker(del_ALT)
    GRanges_delecje <- GRanges(seqnames = delecje$seqnames, IRanges(start = delecje$start, end = delecje$end))
    seqlevelsStyle(GRanges_delecje) <- "UCSC"
    seqlevelsStyle(genome) <- "UCSC"

    # Analyze sequences before and after deletions
    deletion_before_after = deletions_before_after(GRanges_delecje, k)
    colnames(deletion_before_after) <- c("seqnames", "start", "end", "width", "class", "before", "indels", "after")
    deletion_before_after$class <- "DEL"

  }

  # Process insertions (INS) if there are any
  if (!is.null(ins_ALT) && nrow(ins_ALT) != 0){
    # Process insertions (INS) if there are any
    insercje <- insertions_maker(ins_ALT)
    GRanges_insercje <- GRanges(seqnames = insercje$seqnames, IRanges(start = insercje$start, end = insercje$end))
    seqlevelsStyle(GRanges_insercje) <- "UCSC"
    seqlevelsStyle(genome) <- "UCSC"

    # Analyze sequences before and after insertions
    insertion_before_after = insertions_before_after(GRanges_insercje, k)
    colnames(insertion_before_after) <- c("seqnames", "start", "end", "width", "class", "before", "indels", "after")
    insertion_before_after$class <- "INS"
  }

  # Combine deletion and insertion data
  if (!is.null(ins_ALT) && !is.null(del_ALT)) {
    indel_before_after <- rbind(deletion_before_after, insertion_before_after)
  }
  if (!is.null(ins_ALT) && is.null(del_ALT)) {
    indel_before_after <- insertion_before_after
  }
  if (!is.null(del_ALT) && is.null(ins_ALT)){
    indel_before_after <- deletion_before_after
  }

  # Count repeats in the sequences
  repeats <- REPEATS_counter(indel_before_after)

  # Count matching nucleotides and classify
  nucleotide <- NUCLEOTIDE_counter(repeats)

  # Classify sequences based on repeat count and matching nucleotides
  Classification_vcf <- classification(nucleotide, N)

  MNV_ReportTable<-data.frame()
  MNV_variant_caller<-data.frame()
  data_Template_switch<-data.frame()

  List<-list(class=Classification_vcf, MNV_MH_dependent=MNV_ReportTable, MNV_variant_caller_detected=MNV_variant_caller, Template_switch=data_Template_switch)

  # Process mnv (MNV if there are any)
  if (!is.null(mnv_ALT) && nrow(mnv_ALT) != 0){

    GRanges_mnv <- GRanges(seqnames = mnv_ALT$seqnames, IRanges(start = mnv_ALT$start, end = mnv_ALT$end))
    seqlevelsStyle(GRanges_mnv) <- "UCSC"
    seqlevelsStyle(genome) <- "UCSC"

    # Analyze sequences before and after MVN
    # This code defines a function called 'mnv_before_after' that extracts sequences before and after insertions in genomic data.
    mnv_before_after = mnv_before_after(GRanges_mnv, k)
    colnames(mnv_before_after) <- c("seqnames", "start", "end", "width", "class", "before", "indels", "after")
    mnv_before_after$class <- "MNV"
    mnv_before_after$REF <- mnv_ALT$REF
    specified_range <- GRanges(
      seqnames = mnv_ALT$seqnames,
      IRanges(start = mnv_ALT$start - k, end = mnv_ALT$end + k),
      indels = mnv_ALT$ALT)
    DSB_create<-create_DSB(specified_range)
    ReportTable_data<-simulate_template_switch(DSB_create)
    ReportTable_data_2<-FilterFunction(ReportTable_data, mnv_before_after, k)

    MNV_variant_caller=mnv_before_after

    List<-list(class=Classification_vcf,MNV_MH_dependent=MNV_ReportTable, MNV_variant_caller_detected=MNV_variant_caller, Template_switch=data_Template_switch)

    ReportTable<-ReportTable_data_2[[1]]

    MNV_ReportTable<-unique(ReportTable_data_2[[2]])

    if (nrow(ReportTable_data_2[[1]]) != 0){

      data_Template_switch<-ReportTable
      MNV_variant_caller=mnv_before_after
      MNV_ReportTable<-MNV_ReportTable
      List<-list(class=Classification_vcf,MNV_MH_dependent=MNV_ReportTable, MNV_variant_caller_detected=MNV_variant_caller, Template_switch=data_Template_switch)

    }

  }


  return(List)


}
