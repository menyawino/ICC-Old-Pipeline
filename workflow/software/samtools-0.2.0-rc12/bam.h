/*  bam.h -- BAM API.

    Copyright (C) 2008-2014 Genome Research Ltd.
    Portions copyright (C) 2010-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef BAM_BAM_H
#define BAM_BAM_H

/*!
  @header

  BAM library provides I/O and various operations on manipulating files
  in the BAM (Binary Alignment/Mapping) or SAM (Sequence Alignment/Map)
  format. It now supports importing from or exporting to SAM, sorting,
  merging, generating pileup, and quickly retrieval of reads overlapped
  with a specified region.

  @copyright Genome Research Ltd.
 */

#ifndef VERSION
#define BAM_VERSION "0.2.0+"
#else
#define BAM_VERSION VERSION
#endif

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"

/*! @abstract BAM file handler */
typedef BGZF *bamFile;
#define bam_open(fn, mode) bgzf_open(fn, mode)
#define bam_dopen(fd, mode) bgzf_fdopen(fd, mode)
#define bam_close(fp) bgzf_close(fp)
#define bam_tell(fp) bgzf_tell(fp)
#define bam_seek(fp, pos, dir) bgzf_seek(fp, pos, dir)

/*! @typedef
  @abstract Structure for the alignment header.
  @field n_targets   number of reference sequences
  @field target_name names of the reference sequences
  @field target_len  lengths of the referene sequences
  @field dict        header dictionary
  @field hash        hash table for fast name lookup
  @field rg2lib      hash table for @RG-ID -> LB lookup
  @field l_text      length of the plain text in the header
  @field text        plain text

  @discussion Field hash points to null by default. It is a private
  member.
 */
typedef bam_hdr_t bam_header_t;

// TODO This flag-formatting functionality does not currently exist in htslib
#define BAM_OFDEC          0
#define BAM_OFHEX          1
#define BAM_OFSTR          2

/*! @abstract defautl mask for pileup */
#define BAM_DEF_MASK (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

/*! @typedef
  @abstract Structure for core alignment information.
  @field  tid     chromosome ID, defined by bam_header_t
  @field  pos     0-based leftmost coordinate
  @field  bin     bin calculated by bam_reg2bin()
  @field  qual    mapping quality
  @field  l_qname length of the query name
  @field  flag    bitwise flag
  @field  n_cigar number of CIGAR operations
  @field  l_qseq  length of the query sequence (read)
 */
// typedef struct { ... } bam1_core_t;

/*! @typedef
  @abstract Structure for one alignment.
  @field  core       core information about the alignment
  @field  l_aux      length of auxiliary data
  @field  data_len   current length of bam1_t::data
  @field  m_data     maximum length of bam1_t::data
  @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux

  @discussion Notes:

   1. qname is zero tailing and core.l_qname includes the tailing '\0'.
   2. l_qseq is calculated from the total length of an alignment block
      on reading or from CIGAR.
   3. cigar data is encoded 4 bytes per CIGAR operation.
   4. seq is nybble-encoded according to bam_nt16_table.
 */
// typedef struct { ... } bam1_t;
// NOTE htslib version doesn't have l_aux; use bam_get_l_aux(b) instead
#ifndef SAMTOOLS_HTSLIB_SUPPRESS_HACKS
// NOTE htslib also renames data_len to l_data; this macro may help or hinder
#define data_len l_data
#endif

typedef hts_itr_t *bam_iter_t;

#define bam1_strand(b) (bam_is_rev((b)))
#define bam1_mstrand(b) (bam_is_mrev((b)))

/*! @function
  @abstract  Get the CIGAR array
  @param  b  pointer to an alignment
  @return    pointer to the CIGAR array

  @discussion In the CIGAR array, each element is a 32-bit integer. The
  lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
  length of a CIGAR.
 */
#define bam1_cigar(b) (bam_get_cigar((b)))

/*! @function
  @abstract  Get the name of the query
  @param  b  pointer to an alignment
  @return    pointer to the name string, null terminated
 */
#define bam1_qname(b) (bam_get_qname((b)))

/*! @function
  @abstract  Get query sequence
  @param  b  pointer to an alignment
  @return    pointer to sequence

  @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
  8 for T and 15 for N. Two bases are packed in one byte with the base
  at the higher 4 bits having smaller coordinate on the read. It is
  recommended to use bam1_seqi() macro to get the base.
 */
#define bam1_seq(b) (bam_get_seq((b)))

/*! @function
  @abstract  Get query quality
  @param  b  pointer to an alignment
  @return    pointer to quality string
 */
#define bam1_qual(b) (bam_get_qual((b)))

/*! @function
  @abstract  Get a base on read
  @param  s  Query sequence returned by bam1_seq()
  @param  i  The i-th position, 0-based
  @return    4-bit integer representing the base.
 */
#define bam1_seqi(s, i) (bam_seqi((s), (i)))

/*! @function
  @abstract  Get auxiliary data
  @param  b  pointer to an alignment
  @return    pointer to the concatenated auxiliary data
 */
#define bam1_aux(b) (bam_get_aux((b)))

/*!
  @abstract Verbose level between 0 and 3; 0 is supposed to disable all
  debugging information, though this may not have been implemented.
 */
#define bam_verbose hts_verbose

/*! @abstract Table for converting a nucleotide character to the 4-bit encoding. */
#define bam_nt16_table seq_nt16_table

/*! @abstract Table for converting a 4-bit encoded nucleotide to a letter. */
#define bam_nt16_rev_table seq_nt16_str

#ifdef __cplusplus
extern "C" {
#endif

    /*********************
     * Low-level SAM I/O *
     *********************/

    /*! @abstract TAM file handler */
    typedef samFile *tamFile;

    /*!
      @abstract   Open a SAM file for reading, either uncompressed or compressed by gzip/zlib.
      @param  fn  SAM file name
      @return     SAM file handler
     */
    static inline tamFile samtools_sam_open(const char *fn) { return sam_open(fn, "r"); }
    #undef  sam_open
    #define sam_open samtools_sam_open

    /*!
      @abstract   Close a SAM file handler
      @param  fp  SAM file handler
     */
    // void sam_close(tamFile fp);

    /*!
      @abstract      Read one alignment from a SAM file handler
      @param  fp     SAM file handler
      @param  header header information (ordered names of chromosomes)
      @param  b      read alignment; all members in b will be updated
      @return        0 if successful; otherwise negative
     */
    // int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b);

    /*!
      @abstract       Read header information from a TAB-delimited list file.
      @param  fn_list file name for the list
      @return         a pointer to the header structure

      @discussion Each line in this file consists of chromosome name and
      the length of chromosome.
     */
    bam_header_t *sam_header_read2(const char *fn_list);

    /*!
      @abstract       Read header from a SAM file (if present)
      @param  fp      SAM file handler
      @return         pointer to header struct; 0 if no @SQ lines available
     */
    static inline bam_header_t *sam_header_read(tamFile fp) { return sam_hdr_read(fp); }

    // Note the distressing cast -- bam_name2id is not thread-safe
    static inline int32_t bam_get_tid(const bam_header_t *header, const char *seq_name) { return bam_name2id((bam_header_t *)header, seq_name); }


    /*********************
     * Low-level BAM I/O *
     *********************/

    /*!
      @abstract Initialize a header structure.
      @return   the pointer to the header structure
     */
    static inline bam_header_t *bam_header_init(void) { return bam_hdr_init(); }

    /*!
      @abstract        Destroy a header structure.
      @param  header  pointer to the header
     */
    static inline void bam_header_destroy(bam_header_t *header) { bam_hdr_destroy(header); }

    /*!
      @abstract   Read a header structure from BAM.
      @param  fp  BAM file handler, opened by bam_open()
      @return     pointer to the header structure

      @discussion The file position indicator must be placed at the
      beginning of the file. Upon success, the position indicator will
      be set at the start of the first alignment.
     */
    static inline bam_header_t *bam_header_read(bamFile fp) { return bam_hdr_read(fp); }

    /*!
      @abstract      Write a header structure to BAM.
      @param  fp     BAM file handler
      @param  header pointer to the header structure
      @return        always 0 currently
     */
    static inline int bam_header_write(bamFile fp, const bam_header_t *header) { return bam_hdr_write(fp, header); }

    /*!
      @abstract   Read an alignment from BAM.
      @param  fp  BAM file handler
      @param  b   read alignment; all members are updated.
      @return     number of bytes read from the file

      @discussion The file position indicator must be
      placed right before an alignment. Upon success, this function
      will set the position indicator to the start of the next
      alignment. This function is not affected by the machine
      endianness.
     */
    // int bam_read1(bamFile fp, bam1_t *b);

    int bam_remove_B(bam1_t *b);

    /*!
      @abstract   Write an alignment to BAM.
      @param  fp  BAM file handler
      @param  b   alignment to write
      @return     number of bytes written to the file
     */
    // int bam_write1(bamFile fp, const bam1_t *b);

    /*! @function
      @abstract  Initiate a pointer to bam1_t struct
     */
//#define bam_init1()

    /*! @function
      @abstract  Free the memory allocated for an alignment.
      @param  b  pointer to an alignment
     */
//#define bam_destroy1(b)

    /*!
      @abstract       Format a BAM record in the SAM format
      @param  header  pointer to the header structure
      @param  b       alignment to print
      @return         a pointer to the SAM string
     */
    char *bam_format1(const bam_header_t *header, const bam1_t *b);

    /*! @abstract     Formats a BAM record and writes it and \n to stdout */
    void bam_view1(const bam_header_t *header, const bam1_t *b);

    /*!
      @abstract       Check whether a BAM record is plausibly valid
      @param  header  associated header structure, or NULL if unavailable
      @param  b       alignment to validate
      @return         0 if the alignment is invalid; non-zero otherwise

      @discussion  Simple consistency check of some of the fields of the
      alignment record.  If the header is provided, several additional checks
      are made.  Not all fields are checked, so a non-zero result is not a
      guarantee that the record is valid.  However it is usually good enough
      to detect when bam_seek() has been called with a virtual file offset
      that is not the offset of an alignment record.
     */
    int bam_validate1(const bam_header_t *header, const bam1_t *b);

    // TODO Parses headers, so not yet implemented in terms of htslib
    const char *bam_get_library(bam_header_t *header, const bam1_t *b);


    /***************
     * pileup APIs *
     ***************/

    /*! @typedef
      @abstract Structure for one alignment covering the pileup position.
      @field  b      pointer to the alignment
      @field  qpos   position of the read base at the pileup site, 0-based
      @field  indel  indel length; 0 for no indel, positive for ins and negative for del
      @field  is_del 1 iff the base on the padded read is a deletion
      @field  level  the level of the read in the "viewer" mode

      @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
      difference between the two functions is that the former does not
      set bam_pileup1_t::level, while the later does. Level helps the
      implementation of alignment viewers, but calculating this has some
      overhead.
     */
    // typedef struct { ... } bam_pileup1_t;

    // typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

    // typedef struct incomplete *bam_plp_t;

    // bam_plp_t bam_plp_init(bam_plp_auto_f read, void *data);
    // int bam_plp_push(bam_plp_t iter, const bam1_t *b);
    // const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    // const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    // void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
    // void bam_plp_reset(bam_plp_t iter);
    // void bam_plp_destroy(bam_plp_t iter);

    // typedef struct incomplete *bam_mplp_t;

    // bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func,  void **data);
    // void bam_mplp_destroy(bam_mplp_t iter);
    // void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
    // int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);

    /*! @typedef
      @abstract    Type of function to be called by bam_plbuf_push().
      @param  tid  chromosome ID as is defined in the header
      @param  pos  start coordinate of the alignment, 0-based
      @param  n    number of elements in pl array
      @param  pl   array of alignments
      @param  data user provided data
      @discussion  See also bam_plbuf_push(), bam_plbuf_init() and bam_pileup1_t.
     */
    typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);

    typedef struct {
        bam_plp_t iter;
        bam_pileup_f func;
        void *data;
    } bam_plbuf_t;

    void bam_plbuf_reset(bam_plbuf_t *buf);
    bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data);
    void bam_plbuf_destroy(bam_plbuf_t *buf);
    int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf);

    int bam_pileup_file(bamFile fp, int mask, bam_pileup_f func, void *func_data);

    struct __bam_lplbuf_t;
    typedef struct __bam_lplbuf_t bam_lplbuf_t;

    void bam_lplbuf_reset(bam_lplbuf_t *buf);

    /*! @abstract  bam_plbuf_init() equivalent with level calculated. */
    bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data);

    /*! @abstract  bam_plbuf_destroy() equivalent with level calculated. */
    void bam_lplbuf_destroy(bam_lplbuf_t *tv);

    /*! @abstract  bam_plbuf_push() equivalent with level calculated. */
    int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *buf);


    /*********************
     * BAM indexing APIs *
     *********************/

    typedef hts_idx_t bam_index_t;

    /*!
      @abstract   Build index for a BAM file.
      @discussion Index file "fn.bai" will be created.
      @param  fn  name of the BAM file
      @return     always 0 currently
     */
    static inline int samtools_bam_index_build(const char *fn) { return bam_index_build(fn, 0); }
    #undef  bam_index_build
    #define bam_index_build samtools_bam_index_build

    /*!
      @abstract   Load index from file "fn.bai".
      @param  fn  name of the BAM file (NOT the index file)
      @return     pointer to the index structure
     */
    // bam_index_t *bam_index_load(const char *fn);

    /*!
      @abstract    Destroy an index structure.
      @param  idx  pointer to the index structure
     */
    static inline void bam_index_destroy(bam_index_t *idx) { hts_idx_destroy(idx); }

    /*! @typedef
      @abstract      Type of function to be called by bam_fetch().
      @param  b     the alignment
      @param  data  user provided data
     */
    typedef int (*bam_fetch_f)(const bam1_t *b, void *data);

    /*!
      @abstract Retrieve the alignments that are overlapped with the
      specified region.

      @discussion A user defined function will be called for each
      retrieved alignment ordered by its start position.

      @param  fp    BAM file handler
      @param  idx   pointer to the alignment index
      @param  tid   chromosome ID as is defined in the header
      @param  beg   start coordinate, 0-based
      @param  end   end coordinate, 0-based
      @param  data  user provided data (will be transferred to func)
      @param  func  user defined function
     */
    int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);

    static inline bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end) { return bam_itr_queryi(idx, tid, beg, end); }
    static inline int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b) { return iter? hts_itr_next(fp, iter, b, 0) : bam_read1(fp, b); }
    static inline void bam_iter_destroy(bam_iter_t iter) { bam_itr_destroy(iter); }

    /*!
      @abstract       Parse a region in the format: "chr2:100,000-200,000".
      @discussion     bam_header_t::hash will be initialized if empty.
      @param  header  pointer to the header structure
      @param  str     string to be parsed
      @param  ref_id  the returned chromosome ID
      @param  begin   the returned start coordinate
      @param  end     the returned end coordinate
      @return         0 on success; -1 on failure
     */
    int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *begin, int *end);


    /**************************
     * APIs for optional tags *
     **************************/

    /*!
      @abstract       Retrieve data of a tag
      @param  b       pointer to an alignment struct
      @param  tag     two-character tag to be retrieved

      @return  pointer to the type and data. The first character is the
      type that can be 'iIsScCdfAZH'.

      @discussion  Use bam_aux2?() series to convert the returned data to
      the corresponding type.
    */
    // uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);

    // int32_t bam_aux2i(const uint8_t *s);
    // float bam_aux2f(const uint8_t *s);
    #define bam_aux2d(s) (bam_aux2f((s)))
    // char bam_aux2A(const uint8_t *s);
    // char *bam_aux2Z(const uint8_t *s);

    // int bam_aux_del(bam1_t *b, uint8_t *s);
    // void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data);
    static inline uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2]) { return bam_aux_get(b, tag); } // an alias of bam_aux_get()


    /*****************
     * Miscellaneous *
     *****************/

    /*!
      @abstract Calculate the rightmost coordinate of an alignment on the
      reference genome.

      @param  c      pointer to the bam1_core_t structure
      @param  cigar  the corresponding CIGAR array (from bam1_t::cigar)
      @return        the rightmost coordinate, 0-based
    */
    static inline uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar) { return c->pos + (c->n_cigar? bam_cigar2rlen(c->n_cigar, cigar) : 1); }

    /*!
      @abstract      Calculate the length of the query sequence from CIGAR.
      @param  c      pointer to the bam1_core_t structure
      @param  cigar  the corresponding CIGAR array (from bam1_t::cigar)
      @return        length of the query sequence
    */
    static inline int32_t samtools_bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar) { return bam_cigar2qlen(c->n_cigar, cigar); }
    #undef  bam_cigar2qlen
    #define bam_cigar2qlen samtools_bam_cigar2qlen

#ifdef __cplusplus
}
#endif

/*!
  @abstract    Calculate the minimum bin that contains a region [beg,end).
  @param  beg  start of the region, 0-based
  @param  end  end of the region, 0-based
  @return      bin
 */
static inline int bam_reg2bin(uint32_t beg, uint32_t end)
{
    return hts_reg2bin(beg, end, 14, 5);
}

/*!
  @abstract     Copy an alignment
  @param  bdst  destination alignment struct
  @param  bsrc  source alignment struct
  @return       pointer to the destination alignment struct
 */
// bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)

/*!
  @abstract     Duplicate an alignment
  @param  src   source alignment struct
  @return       pointer to the destination alignment struct
 */
// bam1_t *bam_dup1(const bam1_t *src)

#endif
