/*
*   NAME        ibdprep
*
*   SYNOPSIS
*
*       ibdprep pedFile idLen sexLen twinidLen gtypeLen xLinked doIndex
*               doMLEfreq doMMSibs [ #loci mrkFile ] [ famidLen ]
*
*   PURPOSE
*
*       Prepare the input files required for computing IBDs. At present,
*       only the FASTLINK-based pairwise method of Curtis and Sham (1995)
*       is supported. However, a future release will support Blangero's
*       MENDEL-based Monte Carlo method as well.
*
*       Optionally prepare the input files required by the MENDEL program
*       allfreq, which computes maximum likelihood estimates of the marker
*       allele frequencies.
*
*       While the pairwise method is valid only for noninbred pedigrees,
*       other types of pedigree loops can be handled by FASTLINK if the
*       user explicitly identifies the individuals at which the loops
*       can be broken. This program will detect any loops in the input
*       pedigree data, and inbreeding will be flagged if present. If no
*       inbreeding loops occur, the set of individuals from which loop-
*       breakers may be selected is listed.
*
*   ARGUMENTS
*
*       pedFile     pedigree data file
*
*       idLen       width of IDs
*       sexLen      width of sex field
*       twinidLen   width of twin ID field
*       gtypeLen    width of genotype field
*
*       xLinked     if "y", marker loci are X-linked
*
*       doIndex     if "y", assign sequential IDs and write the indexed
*                   pedigree data to the PEDSYS file pedindex.out
*
*       doMLEfreq   if "y", create the input files for program allfreq, a
*                   MENDEL program which computes MLEs of the marker allele
*                   frequencies
*
*       doMMSibs    if "y", create the input files for MAPMAKER/SIBS
*       #loci       number of marker loci (MAPMAKER/SIBS option only)
*       mrkFile     marker location file (MAPMAKER/SIBS option only)
*
*       famidLen    width of family ID (to be used only if individual IDs
*                   are not unique across pedigrees)
*
*   INPUT
*
*       Pedigree information will be read from a file which must contain
*       the following fixed-width fields: ego, father, and mother IDs,
*       sex, twin ID, and marker genotype. Missing values are specified
*       with a blank field.
*
*       IDs must be alphanumeric, and must be unique across the entire
*       data set. Founders will have blanks in the father and mother ID
*       fields. If one parent is missing, the other parent must also be
*       missing. In the case the pedigree data is already indexed, i.e.
*       the IDs are sequential integers as would be assigned by the PEDSYS
*       program INDEX, "n" should be specified for the doIndex command
*       line argument so that the data set will not be re-indexed. In this
*       case, missing parents should have IDs of zero.
*
*       Sex must be coded 1, M, or m for males and 2, F, or f for females.
*       The missing value for sex is 0, U, u, or blank.
*
*       The twin ID field is used to flag genetically identical individuals.
*       Each set of such individuals must be assigned a unique identifier,
*       where the identifiers are sequential integers beginning with one.
*       This identifier must appear in the twin ID field of every individual
*       in the set.
*
*       Marker genotypes must conform to the PEDSYS nomenclature for
*       autosomal codominant markers. For X-linked markers, genotypes for
*       males must be coded as though the male was homozygous at the marker.
*
*   OUTPUT
*
*       The FASTLINK input file datafile.dat is created. It contains estimated
*       marker allele frequencies calculated by counting the alleles present
*       in the pedigree-data file. The file ped.raw is created, which must
*       then be input to the LINKAGE program makeped to create the FASTLINK
*       input file pedfile.dat.
*
*   AUTHOR
*
*       Thomas Dyer
*       Department of Genetics
*       Southwest Foundation for Biomedical Research
*
*       Copyright (C) 1997 Southwest Foundation for Biomedical Research
*       All rights reserved. Absolutely no warranty, express or implied.
*
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#define MIDLEN 36
#define MAXIND 210000	/* cannot exceed 2^15 - 1 */
#define MAXFAM 210000	/* cannot exceed 2^15 - 1 */
#define MAXPED 210000	/* cannot exceed 2^15 - 1 */
#define MXTWIN 210000
#define MAXLOC 3000
#define MMRKNM 20

#define MSXLEN 2
#define MGTLEN 20
#define MAXALL 500

#define MXGENO MAXALL*(MAXALL+1)/2

#define TRUE 1
#define FALSE 0

#define ERRFILE "ibdprep.err"
#define WRNFILE "ibdprep.wrn"

#define max(a, b)       ( (a) > (b) ? (a) : (b) )
#define min(a, b)       ( (a) < (b) ? (a) : (b) )

#define fclose   (void) fclose
#define fprintf  (void) fprintf
#define memset   (void) memset
#define printf   (void) printf
#define sprintf  (void) sprintf
#define strcat   (void) strcat
#define strcpy   (void) strcpy
#define strncpy  (void) strncpy
#define unlink   (void) unlink

int doIndex;
int doMLEfreq = FALSE;
int doMMSibs = FALSE;
int doMCarlo = FALSE;
int doLinkage = FALSE;

int xLinked = FALSE;
int inMrkFile = FALSE;
int dropSingles = FALSE;
int isInbred = FALSE;

struct Ind {
    char id[MIDLEN+1];
    char pid[MIDLEN+1];
    struct Fam *fam;
    struct Ind *sib;
    char sex;
    char twinid[MIDLEN+1];
    int itwinid;
    char hhid[MIDLEN+1];
    int mrkall[MAXLOC][2];
    int ped;
    int gen;
    int seq;
};

struct Fam {
    struct Ind *fa;
    struct Ind *mo;
    struct Ind *kid1;
    struct Fam *next;
    int nkid;
    int ped;
    int seq;
};

struct Ped {
    struct Fam *fam1;
    int nfam;
    int nind;
    int nfou;
    int seq1;
    int inbred;
    int hasloops;
    int nlbrk;
    int lbrkind;
};

struct Link {
    int ind;
    int fam;
    struct Link *next;
};

struct Loc {
    char mrkName[MMRKNM+1];
    int noLocInfo;
    char *allList[MAXALL];
    double allFreq[MAXALL];
    int allSort[MAXALL];
    int allNumeric;
    int numAll;
    int numTyp;
    int numFouTyp;
};

struct Twin {
    char twinid[MIDLEN+1];
    char sex;
    struct Fam *fam;
    int mrkall[2];
};

typedef struct Ind Ind;
typedef struct Fam Fam;
typedef struct Ped Ped;
typedef struct Link Link;
typedef struct Loc Loc;
typedef struct Twin Twin;

Ind *IndArray[MAXIND];
int IndSort[MAXIND];      /* Ind's sorted lexicographically by id  */
int IndSeq[MAXIND];       /* Ind's in sequentially indexed order   */
int PidSort[MAXIND];      /* Ind's sorted lexicographically by pid */
int NumInd;
int NumFou;
int MaxLbrk;

Fam *FamArray[MAXFAM];
int NumFam;

Ped *PedArray[MAXPED];
int NumPed;

Loc *LocArray[MAXLOC];
int NumLoc;

Twin *TwinArray[MXTWIN];
int NumTwin;

int FamidLen = 0;
int IdLen;
int SexLen;
int TwinIdLen;
int TwinOutLen = 3;
int HHIdLen;
int PidLen;
int GtypeLen;

char PedFile[1024];
char MrkFile[1024];
char LocFile[1024];
char MapFile[1024];

FILE *WrnFP;
char WrnMsg[1024];
int WrnCnt;

FILE *ErrFP;
char ErrMsg[1024];
int ErrCnt;

void addLink (Link**, int*, int*, int, int, int);
void *allocMem (size_t);
void assignSeq (void);
void calcKin2 (void);
void checkLooping (void);
void checkTwins (void);
void cntAlleles (int, char**, int*, int*);
void displayUsage (void);
void fatalError (void);
int findAllele (int, char*);
int findBreaks(int*, int*, char*);
int findInd (char*);
int findPid (char*);
int findTwin (char*);
int getAlleles (char*, char**, int*);
void getCmdLine (int, char**);
void getLocInfo (void);
void getMrkData (void);
void getPedData (int*, char**, int*);
void logWarning (void);
void logError (void);
void makeDir (char*, mode_t);
int makeFams (char**, int*, int);
void makeHHoldMat (void);
void makeLinks (int, Link**, int*, int*);
void makePeds (void);
FILE *openFile (char*, char*);
void point (int, int, int*, int*);
void qSort (char**, int, int, int*, int);
void rmLink(Link**, int*, int*, int, int);
int sameGtype (int*, int*);
void sortInds (void);
void sortPids (void);
void trace (int**, int, int, int*, int*);
int unknown (char*);
void warshall (unsigned char**, int);
void writeLocInfo (void);
void writeMLEfreqFiles (int);
void writeMCarloFiles (int);
void writeIndex (void);
void writeLinkageFiles (int);
void writeMakepedCmd (int);
void writeMMSibsFiles (void);
void writeInfo (void);

int main (int argc, char **argv)
{
    int i, nfam;
    int idFam[MAXIND];
    char *famList[MAXIND];

    getCmdLine(argc, argv);

    WrnFP = openFile(WRNFILE, "w");
    WrnCnt = 0;

    ErrFP = openFile(ERRFILE, "w");
    ErrCnt = 0;
    
    getPedData(&nfam, famList, idFam);
    sortInds();
    if (makeFams(famList, idFam, nfam)) {
        sortInds();
        makeFams(famList, idFam, nfam);
    }
    for (i = 0; i < NumFam; i++)
        free(famList[i]);
    checkTwins();
    if (!doIndex) {
        getLocInfo();
        sortPids();
        getMrkData();
        writeLocInfo();
    }
    makePeds();
    checkLooping();     /* must be done before assignSeq!! */
    assignSeq();
    if (doIndex) {
        calcKin2();		/* inbreeding detected here */
        if (HHIdLen)
            makeHHoldMat();
        writeIndex();
    }

    if (doMLEfreq) {
        for (i = 0; i < NumLoc; i++) {
            writeMLEfreqFiles(i);
        }
    }

    if (doMMSibs)
        writeMMSibsFiles();
    else if (doLinkage) {
        for (i = 0; i < NumLoc; i++) {
            writeLinkageFiles(i);
            if (MaxLbrk <= 1)
                writeMakepedCmd(i);
            writeMCarloFiles(i);	/* set up files for both methods */
        }
    }
    else if (doMCarlo) {
        for (i = 0; i < NumLoc; i++) {
            writeMCarloFiles(i);
        }
    }

    writeInfo();

    if (WrnCnt) {
        fprintf(stdout, "%d warnings were written to file \"%s\".",
                WrnCnt, WRNFILE);
    } else {
        unlink(WRNFILE);
    }

    if (!ErrCnt)
        unlink(ERRFILE);

    return 0;
}

void getCmdLine (int argc, char **argv)
{
    FILE *infp;

    if (argc < 7 || argc > 12)
        displayUsage();

    if (argv[1][0] == 'y' || argv[1][0] == 'Y') {
        doIndex = TRUE;
    }
    else if (argv[1][0] == 'n' || argv[1][0] == 'N') {
        doIndex = FALSE;
        doMLEfreq = TRUE;
    }
    else {
        sprintf(ErrMsg, "doIndex? must be y or n");
        fatalError();
    }
 
    if (doIndex) {
        if (argc > 8)
            displayUsage();

        infp = fopen(argv[2], "r");
        if (!infp) {
            sprintf(ErrMsg, "cannot open pedigree-data file \"%s\"", argv[2]);
            fatalError();
        }
        strcpy(PedFile, argv[2]);
        fclose(infp);
 
        if (sscanf(argv[3], "%d", &IdLen) != 1 || IdLen <= 0) {
            sprintf(ErrMsg, "invalid idLen \"%s\"", argv[3]);
            fatalError();
        }
 
        if (IdLen > MIDLEN) {
            sprintf(ErrMsg, "idLen too large, MIDLEN = %d", MIDLEN);
            fatalError();
        }
 
        if (sscanf(argv[4], "%d", &SexLen) != 1 || SexLen <= 0) {
            sprintf(ErrMsg, "invalid sexLen \"%s\"", argv[4]);
            fatalError();
        }
 
        if (SexLen > MSXLEN) {
            sprintf(ErrMsg, "sexLen too large, MSXLEN = %d", MSXLEN);
            fatalError();
        }
 
        if (sscanf(argv[5], "%d", &TwinIdLen) != 1 || TwinIdLen < 0) {
            sprintf(ErrMsg, "invalid twinidLen \"%s\"", argv[5]);
            fatalError();
        }
 
        if (TwinIdLen > MIDLEN) {
            sprintf(ErrMsg, "twinidLen too large, MIDLEN = %d", MIDLEN);
            fatalError();
        }
 
        if (sscanf(argv[6], "%d", &HHIdLen) != 1 || HHIdLen < 0) {
            sprintf(ErrMsg, "invalid hhidLen \"%s\"", argv[6]);
            fatalError();
        }
 
        if (HHIdLen > MIDLEN) {
            sprintf(ErrMsg, "hhidLen too large, MIDLEN = %d", MIDLEN);
            fatalError();
        }

        if (argc == 8) {
            if (sscanf(argv[7], "%d", &FamidLen) != 1 || FamidLen <= 0) {
                sprintf(ErrMsg, "invalid famidLen \"%s\"", argv[7]);
                fatalError();
            }
            if (FamidLen + IdLen > MIDLEN) {
                sprintf(ErrMsg, "famidLen+idLen too large, MIDLEN = %d",
                        MIDLEN);
                fatalError();
            }
        }
    }
 
    else {
        if (argc < 10)
            displayUsage();

        infp = fopen("pedindex.out", "r");
        if (!infp) {
            sprintf(ErrMsg, "cannot open indexed-pedigree file \"%s\"",
                    "pedindex.out");
            fatalError();
        }
        strcpy(PedFile, "pedindex.out");
        fclose(infp);

        infp = fopen(argv[2], "r");
        if (!infp) {
            sprintf(ErrMsg, "cannot open marker-data file \"%s\"", argv[2]);
            fatalError();
        }
        strcpy(MrkFile, argv[2]);
        fclose(infp);
 
        if (sscanf(argv[3], "%d", &IdLen) != 1 || IdLen <= 0) {
            sprintf(ErrMsg, "invalid idLen \"%s\"", argv[3]);
            fatalError();
        }
 
        if (IdLen > MIDLEN) {
            sprintf(ErrMsg, "idLen too large, MIDLEN = %d", MIDLEN);
            fatalError();
        }
 
        if (sscanf(argv[4], "%d", &GtypeLen) != 1 || GtypeLen < 0) {
            sprintf(ErrMsg, "invalid gtypeLen \"%s\"", argv[4]);
            fatalError();
        }
 
        if (GtypeLen > MGTLEN) {
            sprintf(ErrMsg, "gtypeLen too large, MGTLEN = %d", MGTLEN);
            fatalError();
        }

        if (argv[5][0] == 'y' || argv[5][0] == 'Y')
            xLinked = TRUE;
        else if (argv[5][0] == 'n' || argv[5][0] == 'N')
            xLinked = FALSE;
        else {
            sprintf(ErrMsg, "xLinked? must be y or n");
            fatalError();
        }

        if (sscanf(argv[6], "%d", &NumLoc) != 1 || NumLoc < 0) {
            sprintf(ErrMsg, "invalid #loci \"%s\"", argv[6]);
            fatalError();
        }
 
        if (NumLoc > MAXLOC) {
            sprintf(ErrMsg, "#loci too large, MAXLOC = %d", MAXLOC);
            fatalError();
        }
 
        infp = fopen(argv[7], "r");
        if (!infp) {
            infp = fopen(argv[7], "w");
            if (!infp) {
                sprintf(ErrMsg, "cannot open locus-info file \"%s\"", argv[7]);
                fatalError();
            }
        }
        strcpy(LocFile, argv[7]);
        fclose(infp);

        if (argv[8][0] == 'y' || argv[8][0] == 'Y')
            doMCarlo = TRUE;
        else if (argv[8][0] == 'n' || argv[8][0] == 'N')
            doMCarlo = FALSE;
        else {
            sprintf(ErrMsg, "doMCarlo? must be y or n");
            fatalError();
        }
        doLinkage = !doMCarlo;

        if (argv[9][0] == 'y' || argv[9][0] == 'Y')
            doMMSibs = TRUE;
        else if (argv[9][0] == 'n' || argv[9][0] == 'N')
            doMMSibs = FALSE;
        else {
            sprintf(ErrMsg, "doMMSibs? must be y or n");
            fatalError();
        }

        if (doMMSibs) {
            infp = fopen(argv[10], "r");
            if (!infp) {
                sprintf(ErrMsg, "cannot open map-data file \"%s\"", argv[10]);
                fatalError();
            }
            strcpy(MapFile, argv[10]);
            fclose(infp);
        }

        if (doMMSibs && argc >= 12) {
            if (sscanf(argv[11], "%d", &FamidLen) != 1 || FamidLen < 0) {
                sprintf(ErrMsg, "invalid famidLen \"%s\"", argv[11]);
                fatalError();
            }
            if (FamidLen + IdLen > MIDLEN) {
                sprintf(ErrMsg, "famidLen+idLen too large, MIDLEN = %d",
                        MIDLEN);
                fatalError();
            }
            if (FamidLen > 0 && argc == 12)
                displayUsage();
            if (argc == 13) {
                if (argv[12][0] == 'y' || argv[12][0] == 'Y')
                    inMrkFile = TRUE;
                else if (argv[12][0] == 'n' || argv[12][0] == 'N')
                    inMrkFile = FALSE;
                else {
                    sprintf(ErrMsg, "inMrkFile? must be y or n");
                    fatalError();
                }
            }
        }
        else if (!doMMSibs && argc >= 11) {
            if (sscanf(argv[10], "%d", &FamidLen) != 1 || FamidLen < 0) {
                sprintf(ErrMsg, "invalid famidLen \"%s\"", argv[10]);
                fatalError();
            }
            if (FamidLen + IdLen > MIDLEN) {
                sprintf(ErrMsg, "famidLen+idLen too large, MIDLEN = %d",
                        MIDLEN);
                fatalError();
            }
            if (FamidLen > 0 && argc == 11)
                displayUsage();
            if (argc == 12) {
                if (argv[11][0] == 'y' || argv[11][0] == 'Y')
                    inMrkFile = TRUE;
                else if (argv[11][0] == 'n' || argv[11][0] == 'N')
                    inMrkFile = FALSE;
                else {
                    sprintf(ErrMsg, "inMrkFile? must be y or n");
                    fatalError();
                }
            }
        }

        PidLen = IdLen;
        if (inMrkFile)
            PidLen += FamidLen;
        SexLen = 1;
        TwinIdLen = 3;
    }
}

void displayUsage (void)
{
    printf("Usage: ibdprep doIndex? ...\n");
    printf("\n  if doIndex? = y\n");
    printf("     ibdprep y pedFile idLen sexLen twinidLen hhidLen [ famidLen ]\n");
    printf("\n  if doIndex? = n\n");
    printf("     ibdprep n mrkFile idLen gtypeLen xLinked? #loci locFile\n");
    printf("             doMCarlo? doMMSibs? [ mapFile ] [ famidLen [ inMrkFile? ] ]\n");
    exit(1);
}

void getPedData (int *nfam, char **famList, int *idFam)
{
    int i, recLen, seqid, fseqid, mseqid, ch1;
    char pedFmt[1024], pedRec[1024];
    char id[MIDLEN+1], fa[MIDLEN+1], mo[MIDLEN+1];
    char pid[MIDLEN+1], junk[15];
    char famid[MIDLEN+1], tid[MIDLEN+1];
    char sex[MSXLEN+1];
    char twinid[MIDLEN+1], hhid[MIDLEN+1];
    char prtid[MIDLEN+15];
    FILE *pedfp;
 
    pedfp = openFile(PedFile, "r");
    if (doIndex) {
        if (FamidLen) {
            if (TwinIdLen && HHIdLen)
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc%%%dc%%%dc",
                        FamidLen, IdLen, IdLen, IdLen, SexLen, TwinIdLen,
                        HHIdLen);
            else if (TwinIdLen)
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc%%%dc", FamidLen,
                        IdLen, IdLen, IdLen, SexLen, TwinIdLen);
            else if (HHIdLen)
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc%%%dc", FamidLen,
                        IdLen, IdLen, IdLen, SexLen, HHIdLen);
            else
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc", FamidLen, IdLen,
                        IdLen, IdLen, SexLen);
        }
        else {
            if (TwinIdLen && HHIdLen)
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc%%%dc", IdLen,
                        IdLen, IdLen, SexLen, TwinIdLen, HHIdLen);
            else if (TwinIdLen)
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc", IdLen, IdLen,
                        IdLen, SexLen, TwinIdLen);
            else if (HHIdLen)
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc%%%dc", IdLen, IdLen,
                        IdLen, SexLen, HHIdLen);
            else
                sprintf(pedFmt, "%%%dc%%%dc%%%dc%%%dc", IdLen, IdLen, IdLen,
                        SexLen);
        }
        recLen = FamidLen + 3*IdLen + SexLen + TwinIdLen + HHIdLen + 1;
    }
    else {
        if (FamidLen)
            sprintf(pedFmt,
                    "%%5d%%1c%%5d%%1c%%5d%%1c%%%dc%%1c%%%dc%%13c%%%dc%%%dc",
                    SexLen, TwinIdLen, FamidLen, IdLen);
        else
            sprintf(pedFmt,
                    "%%5d%%1c%%5d%%1c%%5d%%1c%%%dc%%1c%%%dc%%13c%%%dc",
                    SexLen, TwinIdLen, IdLen);
        recLen = FamidLen + IdLen + SexLen + TwinIdLen + 33;
    }
    NumInd = 0;
    NumFou = 0;

    *nfam = 0;
    while (fgets(pedRec, sizeof(pedRec), pedfp)) {
        if (strlen(pedRec) != recLen) {
            sprintf(ErrMsg,
                    "incorrect record length, line %d of pedigree-data file",
                    NumInd + 1);
            fatalError();
        }

        if (NumInd == MAXIND) {
            sprintf(ErrMsg, "too many individuals, MAXIND = %d", MAXIND);
            fatalError();
        }

        if (doIndex) {
            if (FamidLen) {
                if (TwinIdLen && HHIdLen)
                    (void) sscanf(pedRec, pedFmt, famid, id, fa, mo, sex,
                                  twinid, hhid);
                else if (TwinIdLen)
                    (void) sscanf(pedRec, pedFmt, famid, id, fa, mo, sex,
                                  twinid);
                else if (HHIdLen)
                    (void) sscanf(pedRec, pedFmt, famid, id, fa, mo, sex,
                                  hhid);
                else
                    (void) sscanf(pedRec, pedFmt, famid, id, fa, mo, sex);
            }
            else {
                if (TwinIdLen && HHIdLen)
                    (void) sscanf(pedRec, pedFmt, id, fa, mo, sex, twinid,
                                  hhid);
                else if (TwinIdLen)
                    (void) sscanf(pedRec, pedFmt, id, fa, mo, sex, twinid);
                else if (HHIdLen)
                    (void) sscanf(pedRec, pedFmt, id, fa, mo, sex, hhid);
                else
                    (void) sscanf(pedRec, pedFmt, id, fa, mo, sex);
            }
        }
        else {
            if (FamidLen)
                (void) sscanf(pedRec, pedFmt, &seqid, junk, &fseqid, junk,
                              &mseqid, junk, sex, junk, twinid, junk, famid,
                              id);
            else
                (void) sscanf(pedRec, pedFmt, &seqid, junk, &fseqid, junk,
                              &mseqid, junk, sex, junk, twinid, junk, id);
            if (--seqid != NumInd || --fseqid > NumInd || --mseqid > NumInd) {
                sprintf(ErrMsg, "pedigree-data file not correctly indexed");
                fatalError();
            }
            strcpy(pid, id);
            pid[PidLen] = '\0';
            memset(fa, ' ', IdLen);
            if (fseqid >= 0)
                strcpy(fa, IndArray[fseqid]->pid);
            memset(mo, ' ', IdLen);
            if (mseqid >= 0)
                strcpy(mo, IndArray[mseqid]->pid);
        }

        id[IdLen] = '\0';
        fa[IdLen] = '\0';
        mo[IdLen] = '\0';

        if (FamidLen) {
            famid[FamidLen] = '\0';
            strcpy(tid, famid);
            strcat(tid, id);
            strcpy(id, tid);
            if (!unknown(fa))
                strcpy(tid, famid);
            else {
                memset(tid, ' ', FamidLen);
                tid[FamidLen] = '\0';
            }
            strcat(tid, fa);
            strcpy(fa, tid);
            if (!unknown(mo))
                strcpy(tid, famid);
            else {
                memset(tid, ' ', FamidLen);
                tid[FamidLen] = '\0';
            }
            strcat(tid, mo);
            strcpy(mo, tid);
        }

        IndArray[NumInd] = (Ind *) allocMem(sizeof(Ind));
        strcpy(IndArray[NumInd]->id, id);
        strcpy(IndArray[NumInd]->pid, pid);

        if (FamidLen)
            sprintf(prtid, "FAMID=\"%s\" ID=\"%s\"", famid, &id[FamidLen]);
        else
            sprintf(prtid, "ID=\"%s\"", &id[FamidLen]);

        sex[SexLen] = '\0';
        i = 0;
        while (i < SexLen && sex[i] == ' ') i++;
        if (i == SexLen) i = 0;
        switch (sex[i]) {
            case '1':
            case 'M':
            case 'm':
                IndArray[NumInd]->sex = 1;
                break;
            case '2':
            case 'F':
            case 'f':
                IndArray[NumInd]->sex = 2;
                break;
            case ' ':
            case '0':
            case 'U':
            case 'u':
                IndArray[NumInd]->sex = 0;
                break;
            default:
                sprintf(ErrMsg,
    "sex must be coded (1,2,0), (M,F,U), or (m,f,u)\n       %s SEX=\"%s\"",
                        prtid, sex);
                logError();
        }

        IndArray[NumInd]->twinid[0] = '\0';
        if (TwinIdLen) {
            twinid[TwinIdLen] = '\0';
            for (i = 0; i < TwinIdLen; i++) {
                if (twinid[i] != ' ' && twinid[i] != '\t' &&
                    twinid[i] != '0') break;
            }
            if (i < TwinIdLen)
                strcpy(IndArray[NumInd]->twinid, twinid);
        }

        IndArray[NumInd]->hhid[0] = '\0';
        if (HHIdLen) {
            hhid[HHIdLen] = '\0';
            for (i = 0; i < HHIdLen; i++) {
                if (hhid[i] != ' ' && hhid[i] != '\t' &&
                    hhid[i] != '0') break;
            }
            if (i < HHIdLen)
                strcpy(IndArray[NumInd]->hhid, hhid);
        }

        if (!unknown(fa) || !unknown(mo)) {
            if (unknown(fa) || unknown(mo)) {
                sprintf(ErrMsg,
    "both parents must be known or unknown\n       %s FA=\"%s\" MO=\"%s\"",
                        prtid, &fa[FamidLen], &mo[FamidLen]);
                logError();
            }
            if (!strcmp(id, fa)) {
                sprintf(ErrMsg,
    "individual has same ID as father\n       %s FA=\"%s\" MO=\"%s\"",
                        prtid, &fa[FamidLen], &mo[FamidLen]);
                logError();
            }
            if (!strcmp(id, mo)) {
                sprintf(ErrMsg,
    "individual has same ID as mother\n       %s FA=\"%s\" MO=\"%s\"",
                        prtid, &fa[FamidLen], &mo[FamidLen]);
                logError();
            }
            if (!strcmp(fa, mo)) {
                sprintf(ErrMsg,
    "father has same ID as mother\n       %s FA=\"%s\" MO=\"%s\"",
                        prtid, &fa[FamidLen], &mo[FamidLen]);
                logError();
            }
            famList[*nfam] = (char *) allocMem((size_t)(2*(IdLen+FamidLen)+1));
            sprintf(famList[*nfam], "%s%s", fa, mo);
            idFam[NumInd] = *nfam;
            IndArray[NumInd]->gen = -1;
            (*nfam)++;
        }
        else {
            idFam[NumInd] = -1;
            IndArray[NumInd]->gen = 0;
            NumFou++;
        }

        IndArray[NumInd]->fam = NULL;
        IndArray[NumInd]->sib = NULL;
        IndArray[NumInd]->ped = -1;
        NumInd++;
    }

    if (ErrCnt) {
        sprintf(ErrMsg, "%d data errors found. See file \"%s\".",
                ErrCnt, ERRFILE);
        fatalError();
    }

    fclose(pedfp);
}

int unknown (char *id)
{
    char *p = id;

    while (*p) {
        if (*p != ' ' && *p != '\t' && *p != '0') return FALSE;
        p++;
    }

    return TRUE;
}

void sortInds (void)
{
    int i;
    char *idList[MAXIND];
    char famid[MIDLEN+1], id[MIDLEN+1], prtid[MIDLEN+15];

    for (i = 0; i < NumInd; i++) {
        idList[i] = (char *) allocMem((size_t)(FamidLen+IdLen+1));
        strcpy(idList[i], IndArray[i]->id);
    }

    qSort(idList, FamidLen+IdLen, NumInd, IndSort, FALSE);

    for (i = 1; i < NumInd; i++) {
        if (!strcmp(IndArray[IndSort[i]]->id,
                    IndArray[IndSort[i-1]]->id))
        {
            if (FamidLen) {
                strcpy(famid, IndArray[IndSort[i]]->id);
                famid[FamidLen] = 0;
                strcpy(id, &IndArray[IndSort[i]]->id[FamidLen]);
                sprintf(prtid, "FAMID=\"%s\" ID=\"%s\"", famid, id);
            }
            else
                sprintf(prtid, "ID=\"%s\"", IndArray[IndSort[i]]->id);
            sprintf(ErrMsg, "individual appears more than once, %s", prtid);
            logError();
        }
    }

    for (i = 0; i < NumInd; i++) {
        IndArray[IndSort[i]]->seq = (int) i;
        free(idList[i]);
    }

    if (ErrCnt) {
        sprintf(ErrMsg, "%d data errors found. See file \"%s\".",
                ErrCnt, ERRFILE);
        fatalError();
    }
}

int makeFams (char **famList, int *idFam, int nfam)
{
    int i, redo = 0;
    int tNumInd;
    int ndx, ord[MAXIND], famndx[MAXIND];
    char fa[MIDLEN+1], mo[MIDLEN+1];
    char famid[MIDLEN+1], id[MIDLEN+1], prtid[MIDLEN+15];
    Ind *indp, *kidp;

    if (!nfam) {
        NumFam = 0;
        return 0;
    }

    qSort(famList, 2*(FamidLen+IdLen), nfam, ord, FALSE);

    /* can't update NumInd until we exit makeFams because
       that would screw up the search in findInd */
    tNumInd = NumInd;

    FamArray[0] = (Fam *) allocMem(sizeof(Fam));
    strncpy(fa, famList[ord[0]], FamidLen+IdLen);
    fa[FamidLen+IdLen] = '\0';

    if (FamidLen) {
        strcpy(famid, fa);
        famid[FamidLen] = 0;
        strcpy(id, &fa[FamidLen]);
        sprintf(prtid, "FAMID=\"%s\" FA=\"%s\"", famid, id);
    }
    else
        sprintf(prtid, "FA=\"%s\"", fa);

    if ((ndx = findInd(fa)) == -1) {
        sprintf(WrnMsg, "record added for father, %s", prtid);
        logWarning();
        IndArray[tNumInd] = (Ind *) allocMem(sizeof(Ind));
        strcpy(IndArray[tNumInd]->id, fa);
        strcpy(IndArray[tNumInd]->pid, fa);
        IndArray[tNumInd]->sex = 1;
        IndArray[tNumInd]->twinid[0] = '\0';
        IndArray[tNumInd]->hhid[0] = '\0';
        IndArray[tNumInd]->gen = 0;
        idFam[tNumInd] = -1;
        NumFou++;
        IndArray[tNumInd]->fam = NULL;
        IndArray[tNumInd]->sib = NULL;
        IndArray[tNumInd]->ped = -1;
        tNumInd++;
        redo = 1;
    }
    else {
        if (IndArray[ndx]->sex != 1) {
            sprintf(WrnMsg, "sex code changed to male for father, %s", prtid);
            logWarning();
            IndArray[ndx]->sex = 1;
        }
        FamArray[0]->fa = IndArray[ndx];
    }

    strncpy(mo, &famList[ord[0]][FamidLen+IdLen], FamidLen+IdLen);
    mo[FamidLen+IdLen] = '\0';

    if (FamidLen) {
        strcpy(famid, mo);
        famid[FamidLen] = 0;
        strcpy(id, &mo[FamidLen]);
        sprintf(prtid, "FAMID=\"%s\" MO=\"%s\"", famid, id);
    }
    else
        sprintf(prtid, "MO=\"%s\"", mo);

    if ((ndx = findInd(mo)) == -1) {
        sprintf(WrnMsg, "record added for mother, %s", prtid);
        logWarning();
        IndArray[tNumInd] = (Ind *) allocMem(sizeof(Ind));
        strcpy(IndArray[tNumInd]->id, mo);
        strcpy(IndArray[tNumInd]->pid, mo);
        IndArray[tNumInd]->sex = 2;
        IndArray[tNumInd]->twinid[0] = '\0';
        IndArray[tNumInd]->hhid[0] = '\0';
        IndArray[tNumInd]->gen = 0;
        idFam[tNumInd] = -1;
        NumFou++;
        IndArray[tNumInd]->fam = NULL;
        IndArray[tNumInd]->sib = NULL;
        IndArray[tNumInd]->ped = -1;
        tNumInd++;
        redo = 1;
    }
    else {
        if (IndArray[ndx]->sex != 2) {
            sprintf(WrnMsg, "sex code changed to female for mother, %s", prtid);
            logWarning();
            IndArray[ndx]->sex = 2;
        }
        FamArray[0]->mo = IndArray[ndx];
    }

    FamArray[0]->kid1 = NULL;
    FamArray[0]->next = NULL;
    FamArray[0]->nkid = 0;
    FamArray[0]->ped = -1;
    famndx[ord[0]] = 0;
    NumFam = 1;

    for (i = 1; i < nfam; i++) {
        if (strcmp(famList[ord[i-1]], famList[ord[i]])) {
            if (NumFam == MAXFAM) {
                sprintf(ErrMsg, "too many families, MAXFAM = %d", MAXFAM);
                fatalError();
            }

            FamArray[NumFam] = (Fam *) allocMem(sizeof(Fam));
            strncpy(fa, famList[ord[i]], FamidLen+IdLen);
            fa[FamidLen+IdLen] = '\0';

            if (FamidLen) {
                strcpy(famid, fa);
                famid[FamidLen] = 0;
                strcpy(id, &fa[FamidLen]);
                sprintf(prtid, "FAMID=\"%s\" FA=\"%s\"", famid, id);
            }
            else
                sprintf(prtid, "FA=\"%s\"", fa);

            if ((ndx = findInd(fa)) == -1) {
                sprintf(WrnMsg, "record added for father, %s", prtid);
                logWarning();
                IndArray[tNumInd] = (Ind *) allocMem(sizeof(Ind));
                strcpy(IndArray[tNumInd]->id, fa);
                strcpy(IndArray[tNumInd]->pid, fa);
                IndArray[tNumInd]->sex = 1;
                IndArray[tNumInd]->twinid[0] = '\0';
                IndArray[tNumInd]->hhid[0] = '\0';
                IndArray[tNumInd]->gen = 0;
                idFam[tNumInd] = -1;
                NumFou++;
                IndArray[tNumInd]->fam = NULL;
                IndArray[tNumInd]->sib = NULL;
                IndArray[tNumInd]->ped = -1;
                tNumInd++;
                redo = 1;
            }
            else {
                if (IndArray[ndx]->sex != 1) {
                    sprintf(WrnMsg,
                            "sex code changed to male for father, %s", prtid);
                    logWarning();
                    IndArray[ndx]->sex = 1;
                }
                FamArray[NumFam]->fa = IndArray[ndx];
            }

            strncpy(mo, &famList[ord[i]][FamidLen+IdLen], FamidLen+IdLen);
            mo[FamidLen+IdLen] = '\0';

            if (FamidLen) {
                strcpy(famid, mo);
                famid[FamidLen] = 0;
                strcpy(id, &mo[FamidLen]);
                sprintf(prtid, "FAMID=\"%s\" MO=\"%s\"", famid, id);
            }
            else
                sprintf(prtid, "MO=\"%s\"", mo);

            if ((ndx = findInd(mo)) == -1) {
                sprintf(WrnMsg, "record added for mother, %s", prtid);
                logWarning();
                IndArray[tNumInd] = (Ind *) allocMem(sizeof(Ind));
                strcpy(IndArray[tNumInd]->id, mo);
                strcpy(IndArray[tNumInd]->pid, mo);
                IndArray[tNumInd]->sex = 2;
                IndArray[tNumInd]->twinid[0] = '\0';
                IndArray[tNumInd]->hhid[0] = '\0';
                IndArray[tNumInd]->gen = 0;
                idFam[tNumInd] = -1;
                NumFou++;
                IndArray[tNumInd]->fam = NULL;
                IndArray[tNumInd]->sib = NULL;
                IndArray[tNumInd]->ped = -1;
                tNumInd++;
                redo = 1;
            }
            else {
                if (IndArray[ndx]->sex != 2) {
                    sprintf(WrnMsg,
                            "sex code changed to female for mother, %s", prtid);
                    logWarning();
                    IndArray[ndx]->sex = 2;
                }
                FamArray[NumFam]->mo = IndArray[ndx];
            }

            FamArray[NumFam]->kid1 = NULL;
            FamArray[NumFam]->next = NULL;
            FamArray[NumFam]->nkid = 0;
            FamArray[NumFam]->ped = -1;
            NumFam++;
        }

        famndx[ord[i]] = NumFam - 1;
    }

    NumInd = tNumInd;

    if (ErrCnt) {
        sprintf(ErrMsg, "%d data errors found. See file \"%s\".",
                ErrCnt, ERRFILE);
        fatalError();
    }

    if (redo) {
        for (i = 0; i < NumFam; i++)
            free(FamArray[i]);
        return 1;
    }

    for (i = 0; i < NumInd; i++) {
        indp = IndArray[i];
        if (idFam[i] >= 0) {
            indp->fam = FamArray[famndx[idFam[i]]];
            if (indp->fam->kid1 == NULL) {
                indp->fam->kid1 = indp;
                indp->fam->nkid++;
            }
            else {
                kidp = indp->fam->kid1;
                while (kidp->sib != NULL)
                    kidp = kidp->sib;
                kidp->sib = indp;
                indp->fam->nkid++;
            }
        }
    }

    return 0;
}

/*
 *  IndArray[ findInd( "id" ) ] = "id"
 */
int findInd (char *id)
{
    int ndx, lo, hi, cmp;

    lo = 0;
    hi = NumInd - 1;
    while (hi >= lo) {
        ndx = (lo + hi)/2;
        cmp = strcmp(id, IndArray[IndSort[ndx]]->id);
        if (cmp == 0)
            return IndSort[ndx];
        else if (cmp < 0)
            hi = ndx - 1;
        else
            lo = ndx + 1;
    }

    return -1;
}

void checkTwins (void)
{
    int i;
    Ind *indp;
    Twin *twinp;

    NumTwin = 0;
    for (i = 0; i < NumInd; i++) {
        indp = IndArray[i];
        if (!strlen(indp->twinid)) {
            indp->itwinid = 0;
            continue;
        }

        indp->itwinid = findTwin(indp->twinid);
        if (indp->itwinid) {
            twinp = TwinArray[indp->itwinid-1];
            if (indp->sex != twinp->sex) {
                sprintf(ErrMsg, "MZ twins of different sex, twin ID = [%s]",
                        indp->twinid);
                logError();
            }
            if (indp->fam != twinp->fam) {
                sprintf(ErrMsg, "MZ twins not in same family, twin ID = [%s]",
                        indp->twinid);
                logError();
            }
        }
        else {
            if (NumTwin == MXTWIN) {
                sprintf(ErrMsg, "too many MZ twins, MXTWIN = %d", MXTWIN);
                logError();
            }

            twinp = (Twin *) allocMem(sizeof(Twin));
            strcpy(twinp->twinid, indp->twinid);
            twinp->sex = indp->sex;
            twinp->fam = indp->fam;

            TwinArray[NumTwin] = twinp;
            NumTwin++;
	    if (NumTwin > 999) {
		TwinOutLen = 5;
	    }

            indp->itwinid = NumTwin;
        }
    }

    if (ErrCnt) {
        sprintf(ErrMsg, "%d data errors found. See file \"%s\".",
                ErrCnt, ERRFILE);
        fatalError();
    }
}

int findTwin (char *twinid)
{
    int i;

    for (i = 0; i < NumTwin; i++) {
        if (!strcmp(twinid, TwinArray[i]->twinid))
            return i + 1;
    }

    return 0;
}

void getLocInfo (void)
{
    int i, loc;
    char locRec[10000], *recp;
    Loc *locp;
    FILE *locfp;

    for (loc = 0; loc < NumLoc; loc++) {
        LocArray[loc] = (Loc *) allocMem(sizeof(Loc));
        locp = LocArray[loc];
        locp->mrkName[0] = '\0';
        locp->noLocInfo = TRUE;
        locp->allNumeric = TRUE;
        locp->numAll = 0;
        locp->numTyp = 0;
        locp->numFouTyp = 0;
    }

    for (i = 0; i < NumInd; i++) {
        for (loc = 0; loc < NumLoc; loc++) {
            IndArray[i]->mrkall[loc][0] = -1;
            IndArray[i]->mrkall[loc][1] = -1;
        }
    }

    locfp = openFile(LocFile, "r");

    loc = 0;
    while (fgets(locRec, sizeof(locRec), locfp)) {
        if (!(recp = strtok(locRec, " \t\n"))) {
            sprintf(ErrMsg, "invalid record, line %d of locus-info file",
                    loc + 1);
            fatalError();
        }

        if (loc >= NumLoc) {
            sprintf(ErrMsg,
                    "too many markers in locus-info file, expected %d",
                    NumLoc);
            fatalError();
        }
        locp = LocArray[loc];

        if (strlen(recp) > MMRKNM) {
            sprintf(ErrMsg, "marker name too long, MMRKNM = %d", MMRKNM);
            fatalError();
        }
        strcpy(locp->mrkName, recp);

        while (recp = strtok(NULL, " \t\n")) {
            if (locp->numAll >= MAXALL) {
                sprintf(ErrMsg,
                        "too many alleles for marker %s, MAXALL = %d",
                        locp->mrkName, MAXALL);
                fatalError();
            }
            if (strlen(recp) > MGTLEN) {
                sprintf(ErrMsg, "allele name too long, MGTLEN = %d", MGTLEN);
                fatalError();
            }

            locp->allList[locp->numAll] =
                         (char *) allocMem((size_t)(GtypeLen+1));
            strcpy(locp->allList[locp->numAll], recp);

            if (!(recp = strtok(NULL, " \t\n")) ||
                    sscanf(recp, "%lf", &locp->allFreq[locp->numAll]) != 1) {
                sprintf(ErrMsg,
                        "invalid record, line %d of locus-info file",
                        loc + 1);
                fatalError();
            }
            locp->numAll++;
        }
        if (locp->numAll > 0)
            locp->noLocInfo = FALSE;
        loc++;
    }

    if (loc < NumLoc) {
        sprintf(ErrMsg,
                "not enough markers in locus-info file, expected %d",
                NumLoc);
        fatalError();
    }

    fclose(locfp);
}

void sortPids (void)
{
    int i;
    char *pidList[MAXIND];

    for (i = 0; i < NumInd; i++) {
        pidList[i] = (char *) allocMem((size_t)(PidLen+1));
        if (inMrkFile)
            strcpy(pidList[i], IndArray[i]->id);
        else
            strcpy(pidList[i], IndArray[i]->pid);
    }

    qSort(pidList, PidLen, NumInd, PidSort, FALSE);
}

void getMrkData (void)
{
    int i, recLen, nrec, seqid, ndx, ch1;
    int loc, gtStart, maxFreq, allCnt[MAXLOC][MAXALL];
    double sumFreq;
    char mrkFmt[100], gtypeFmt[100], mrkRec[100000], cnum[10];
    char pid[MIDLEN+1];
    char famid[MIDLEN+1], id[MIDLEN+1], prtid[MIDLEN+15];
    char gtype[MAXLOC][MGTLEN+1], all1[MGTLEN+1], all2[MGTLEN+1];
    char *allele[2];
    Ind *indp;
    Loc *locp;
    Twin *twinp;
    FILE *mrkfp;
 
    mrkfp = openFile(MrkFile, "r");

    allele[0] = all1;
    allele[1] = all2;

    sprintf(mrkFmt, "%%%dc", PidLen);
    sprintf(gtypeFmt, "%%%dc", GtypeLen);
    recLen = PidLen + NumLoc*GtypeLen + 1;

    nrec = 0;
    while (fgets(mrkRec, sizeof(mrkRec), mrkfp)) {
        nrec++;
        if (strlen(mrkRec) != recLen) {
            sprintf(ErrMsg,
                    "incorrect record length, line %d of marker-data file",
                    nrec);
            fatalError();
        }

        (void) sscanf(mrkRec, mrkFmt, pid);
        pid[PidLen] = '\0';

        gtStart = PidLen;
        for (i = 0; i < NumLoc; i++) {
            (void) sscanf(mrkRec+gtStart, gtypeFmt, gtype[i]);
            gtype[i][GtypeLen] = '\0';
            gtStart += GtypeLen;
        }
        if ((ndx = findPid(pid)) == -1) {
            continue;
        }
        else
            indp = IndArray[ndx];

        if (FamidLen && inMrkFile) {
            strcpy(famid, pid);
            famid[FamidLen] = 0;
            sprintf(prtid, "FAMID=\"%s\" ID=\"%s\"", famid, &pid[FamidLen]);
        }
        else
            sprintf(prtid, "ID=\"%s\"", pid);

        for (loc = 0; loc < NumLoc; loc++) {
            locp = LocArray[loc];
            if (getAlleles(gtype[loc], allele, &locp->allNumeric)) {
                if (xLinked && indp->sex == 2 &&
                    (!strlen(allele[0]) && strlen(allele[1]) ||
                      strlen(allele[0]) && !strlen(allele[1])))
                {
                    sprintf(ErrMsg,
                "invalid female genotype at marker %s\n       %s Gtype=\"%s\"",
                            locp->mrkName, prtid, gtype[loc]);
                    logError();
                }
                else if (xLinked && indp->sex == 1 &&
                         strlen(allele[0]) && strlen(allele[1]) &&
                         strcmp(allele[0], allele[1]))
                {
                    sprintf(ErrMsg,
                "invalid male genotype at marker %s\n       %s Gtype=\"%s\"",
                            locp->mrkName, prtid, gtype[loc]);
                    logError();
                }
                else {
                    if (xLinked && indp->sex == 1) {
                        if (!strlen(allele[0]) && strlen(allele[1]))
                            strcpy(allele[0], allele[1]);
                        else if (strlen(allele[0]) && !strlen(allele[1]))
                            strcpy(allele[1], allele[0]);
                    }
                    cntAlleles(loc, allele, allCnt[loc], indp->mrkall[loc]);
                    if (indp->mrkall[loc][1] != -1) {
                        locp->numTyp++;
                        if (indp->fam == NULL)
                            locp->numFouTyp++;
                    }
                }
            }
            else {
                sprintf(ErrMsg,
                "invalid genotype at marker %s\n       %s Gtype=\"%s\"",
                        locp->mrkName, prtid, gtype[loc]);
                logError();
            }
        }
    }

    for (loc = 0; loc < NumLoc; loc++) {
        locp = LocArray[loc];
        for (i = 0; i < NumTwin; i++) {
            twinp = TwinArray[i];
            twinp->mrkall[0] = -1;
            twinp->mrkall[1] = -1;
        }
        for (i = 0; i < NumInd; i++) {
            indp = IndArray[i];
            if (indp->itwinid) {
                twinp = TwinArray[indp->itwinid-1];
                if (twinp->mrkall[0] != -1) {
                    if (indp->mrkall[loc][0] != -1 &&
                        !sameGtype(indp->mrkall[loc], twinp->mrkall))
                    {
                        sprintf(ErrMsg,
            "MZ twins have different genotypes at marker %s, twin ID = [%s]",
                                locp->mrkName, twinp->twinid);
                        logError();
                    }
                }
                else {
                    twinp->mrkall[0] = indp->mrkall[loc][0];
                    twinp->mrkall[1] = indp->mrkall[loc][1];
                }
            }
        }
    }

    if (ErrCnt) {
        sprintf(ErrMsg, "%d data errors found. See file \"%s\".",
                ErrCnt, ERRFILE);
        fatalError();
    }

    for (loc = 0; loc < NumLoc; loc++) {
        locp = LocArray[loc];

        if (locp->noLocInfo) {
            maxFreq = 0;
            sumFreq = 0.;
            for (i = 0; i < locp->numAll; i++) {
                locp->allFreq[i] = allCnt[loc][i] / (2. * (double)locp->numTyp);
                sprintf(cnum, "%8.6f", locp->allFreq[i]);
                locp->allFreq[i] = atof(cnum);
                if (locp->allFreq[i] > locp->allFreq[maxFreq])
                    maxFreq = i;
                sumFreq += locp->allFreq[i];
            }
            locp->allFreq[maxFreq] = locp->allFreq[maxFreq] - sumFreq + 1;
        }

        qSort(locp->allList, GtypeLen, (int) locp->numAll, locp->allSort,
              locp->allNumeric);
    }

    fclose(mrkfp);
}

int sameGtype (int *mrkall1, int *mrkall2)
{
    if (mrkall1[0] == mrkall2[0] && mrkall1[1] == mrkall2[1])
        return 1;
    if (mrkall1[1] == mrkall2[0] && mrkall1[0] == mrkall2[1])
        return 1;
    return 0;
}

/*
 *  IndArray[ findPid( "pid" ) ] = "pid"
 */
int findPid (char *pid)
{
    int ndx, lo, hi, cmp;

    lo = 0;
    hi = NumInd - 1;
    while (hi >= lo) {
        ndx = (lo + hi)/2;
        if (inMrkFile)
            cmp = strcmp(pid, IndArray[PidSort[ndx]]->id);
        else
            cmp = strcmp(pid, IndArray[PidSort[ndx]]->pid);
        if (cmp == 0)
            return PidSort[ndx];
        else if (cmp < 0)
            hi = ndx - 1;
        else
            lo = ndx + 1;
    }

    return -1;
}

void writeLocInfo (void)
{
    int i, loc;
    Loc *locp;
    FILE *locfp;
 
    locfp = openFile(LocFile, "w");

    for (loc = 0; loc < NumLoc; loc++) {
        locp = LocArray[loc];
        if (!strlen(locp->mrkName))
            sprintf(locp->mrkName, "marker%d", loc + 1);
        fprintf(locfp, "%s", locp->mrkName);
        for (i = 0; i < locp->numAll; i++)
            fprintf(locfp, " %s %8.6f", locp->allList[locp->allSort[i]],
                    locp->allFreq[locp->allSort[i]]);
        fprintf(locfp, "\n");
    }

    fclose(locfp);
}

void makePeds (void)
{
    int i, j, curind, ind;
    int k, m, n, ok, ip, isave, *perm;
    char famid[MIDLEN+1], prtid[MIDLEN+15];
    int genfnd, lastgen, nped, fgen, mgen;
    int *relate[5], *stack, *state;
    Fam **fampp;
    Ind *indp;

    perm = (int *) allocMem((size_t)(NumInd*sizeof(int)));
    m = 0;
    n = NumInd - 1;
    for (i = 0; i < NumInd; i++) {
        indp = IndArray[IndSort[i]];
        if (indp->fam) {
            perm[n] = i + 2*NumInd;
/*printf("n=%d perm[n]=%d\n", n, perm[n]);fflush(stdout);*/
            n--;
        }
        else {
            perm[m] = i;
/*printf("m=%d perm[m]=%d\n", m, perm[m]);fflush(stdout);*/
            m++;
        }
    }

    m = 0;
    while (1) {
        ok = 0;
        for (k = m; k < NumInd; k++) {
            if (perm[k] < NumInd) {
                ok = 1;
                break;
            }
        }
        if (!ok) {
            ip = perm[m] % NumInd;
/*printf("ip=%d m=%d perm[m]=%d id=%s\n", ip, m, perm[m], IndArray[IndSort[ip]]->id);fflush(stdout);*/
            if (FamidLen) {
                strcpy(famid, IndArray[IndSort[ip]]->id);
                famid[FamidLen] = 0;
                sprintf(prtid, "FAMID=\"%s\" ID=\"%s\"", famid, &IndArray[IndSort[ip]]->id[FamidLen]);
            }
            else
                sprintf(prtid, "ID=\"%s\"", IndArray[IndSort[ip]]->id);
            sprintf(ErrMsg, "an individual near %s is his/her own ancestor", prtid);
            fatalError();
        }
        isave = perm[k];
/*printf("isave=%d m=%d perm[m]=%d k=%d perm[k]=%d\n", isave, m, perm[m], k, perm[k]);fflush(stdout);*/
        perm[k] = perm[m];
        perm[m] = isave;
        m++;
        if (m == NumInd - 1)
            break;
        for (i = m; i < NumInd; i++) {
            ip = perm[i] % NumInd;
            indp = IndArray[IndSort[ip]];
            if (indp->fam && indp->fam->fa == IndArray[IndSort[isave]])
                perm[i] -= NumInd;
            if (indp->fam && indp->fam->mo == IndArray[IndSort[isave]])
                perm[i] -= NumInd;
/*printf("m=%d i=%d perm[i]=%d ip=%d\n", m, i, perm[i], ip);fflush(stdout);*/
        }
    }

    free(perm);

    genfnd = NumFou;
    while (genfnd < NumInd) {
        lastgen = genfnd;
        for (i = 0; i < NumInd; i++) {
            indp = IndArray[IndSort[i]];
            if (indp->gen < 0) {
                if (indp->fam) {
                    fgen = indp->fam->fa->gen;
                    mgen = indp->fam->mo->gen;
                }
                else {
                    fgen = 0;
                    mgen = 0;
                }
                if (fgen >= 0 && mgen >= 0) {
                    indp->gen = max(fgen,mgen) + 1;
                    genfnd++;
                }
            }
        }
        if (genfnd == lastgen) {
            sprintf(ErrMsg, "pedigree error detected while assigning generation numbers");
            fatalError();
        }
    }

    for (i = 0; i < 5; i++) {
        relate[i] = (int *) allocMem((size_t)(NumInd*sizeof(int)));
        for (ind = 0; ind < NumInd; ind++)
            relate[i][ind] = -1;
    }

    for (curind = 0; curind < NumInd; curind++) {
        indp = IndArray[curind];
        if (indp->fam) {
            ind = findInd(indp->fam->fa->id);
            relate[0][curind] = ind;
            point(curind, ind, relate[2], relate[4]);
            ind = findInd(indp->fam->mo->id);
            relate[1][curind] = ind;
            point(curind, ind, relate[3], relate[4]);
        }
        else {
            relate[0][curind] = -1;
            relate[1][curind] = -1;
        }
    }

    stack = (int *) allocMem((size_t)(NumInd*sizeof(int)));
    state = (int *) allocMem((size_t)(NumInd*sizeof(int)));
    for (ind = 0; ind < NumInd; ind++)
         state[ind] = -1;
    for (ind = 0; ind < NumInd; ind++)
         IndArray[ind]->ped = -1;

    nped = 0;
    for (curind = 0; curind < NumInd; curind++) {
        if (relate[0][curind] == -1 && relate[1][curind] == -1
                && relate[2][curind] == -1 && relate[3][curind] == -1
                && relate[4][curind] == -1)
            continue;
        if (IndArray[curind]->ped == -1) {
            trace(relate, curind, nped, stack, state);
            nped++;
        }
    }

    for (i = 0; i < 5; i++)
        free(relate[i]);
    free(stack);
    free(state);

    for (i = 0; i < nped; i++) {
        PedArray[i] = (Ped *) allocMem(sizeof(Ped));
        PedArray[i]->nfam = 0;
        PedArray[i]->nind = 0;
        PedArray[i]->nfou = 0;
        fampp = &(PedArray[i]->fam1);
        for (j = 0; j < NumFam; j++) {
            if (FamArray[j]->ped == (int) i) {
                FamArray[j]->seq = PedArray[i]->nfam;
                PedArray[i]->nfam++;
                *fampp = FamArray[j];
                fampp = &(FamArray[j]->next);
            }
        }
    }

    NumPed = nped;
    for (i = 0; i < NumInd; i++) {
        indp = IndArray[i];
        if (indp->ped == -1) {
/*            printf("WARNING: Individual \"%s\" is unrelated.\n", indp->id);*/
            PedArray[NumPed] = (Ped *) allocMem(sizeof(Ped));
            PedArray[NumPed]->nind = 1;
            PedArray[NumPed]->nfam = 0;
            PedArray[NumPed]->nfou = 1;
            PedArray[NumPed]->fam1 = NULL;
            indp->ped = NumPed;
            NumPed++;
        }
        else {
            PedArray[indp->ped]->nind++;
            if (!indp->fam)
                PedArray[indp->ped]->nfou++;
        }
    }
}

void trace (int **relate, int curind, int curped, int *stack, int *state)
{
    int pstack;

    pstack = -1;
    for ( ; ; ) {
        if (curind == -1 || IndArray[curind]->ped != -1)
            curind = stack[pstack];
        else {
            pstack++;
            stack[pstack] = curind;
            IndArray[curind]->ped = curped;
            if (IndArray[curind]->fam)
                IndArray[curind]->fam->ped = curped;
        }
        state[curind]++;
        if (state[curind] > 5) {
            pstack--;
            if (pstack < 0) return;
        }
        else if (state[curind] < 5)
            curind = relate[state[curind]][curind];
        else if (state[curind] == 5 && relate[4][curind] != -1) {
            if (relate[0][relate[4][curind]] == curind)
                curind = relate[1][relate[4][curind]];
            else
                curind = relate[0][relate[4][curind]];
        }
    }
}

void point (int curind, int ind, int *r_next, int *r_off1)
{
    if (r_off1[ind] == -1)
        r_off1[ind] = curind;
    else {
        ind = r_off1[ind];
        while (r_next[ind] != -1)
            ind = r_next[ind];
        r_next[ind] = curind;
    }
}

void checkLooping (void)
{
    int i, j, narcs, maxlbrk;
    Ped *pedp;
    Fam *famp;

    int nlink[MAXFAM];
    int linkInd[MAXIND];
    char lbrkId[MIDLEN+1];
    Link *linkList[MAXFAM];
    Link *linkp, *lastp;

    MaxLbrk = 0;
    for (i = 0; i < NumPed; i++) {
        pedp = PedArray[i];

        narcs = 0;
        famp = pedp->fam1;
        while (famp) {
            narcs += famp->nkid + 2;
            famp = famp->next;
        }

        if (narcs < pedp->nind + pedp->nfam) {
            pedp->hasloops = FALSE;
            pedp->nlbrk = 0;
            continue;
        }

/*        printf("Pedigree %d contains 1 or more loops: ", i+1);*/
        pedp->hasloops = TRUE;

        makeLinks(i, linkList, nlink, linkInd);
        pedp->nlbrk = findBreaks(nlink, linkInd, lbrkId);
        MaxLbrk = max(pedp->nlbrk, MaxLbrk);

        if (pedp->nlbrk == 1) {
            for (j = 0; j < NumInd; j++) {
                if (linkInd[j] && IndArray[IndSort[j]]->fam
                        && nlink[IndArray[IndSort[j]]->fam->seq]) {
                    pedp->lbrkind = IndSort[j];
                    break;
                }
            }
        }
/*
        printf("%d loop-breaker(s) required.\n", pedp->nlbrk);
        if (pedp->nlbrk >= 1) {
            printf("The following individuals are possible loop-breakers:\n");
            for (j = 0; j < NumInd; j++) {
                if (linkInd[j] && IndArray[IndSort[j]]->fam
                        && nlink[IndArray[IndSort[j]]->fam->seq]) {
                    printf("\"%s\"\n", IndArray[IndSort[j]]->id);
                }
            }
        }
*/
        for (j = 0; j < NumFam; j++) {
            linkp = linkList[j];
            while (linkp) {
                lastp = linkp;
                linkp = linkp->next;
                free(lastp);
            }
        }
    }
}

void makeLinks (int ped, Link **linkList, int *nlink, int *linkInd)
{
    int i, j, done;
    Ped *pedp;
    Fam *famp, *famp2;
    Link *linkp, *lastp;

    pedp = PedArray[ped];
    for (i = 0; i < NumFam; i++) {
        linkList[i] = NULL;
        nlink[i] = 0;
    }

    for (i = 0; i < NumInd; i++)
        linkInd[i] = 0;

/*printf("\n");*/
    famp = pedp->fam1;
    while (famp) {
/*printf("fam = %d  fa = %s  fseq = %d  mo = %s  mseq = %d\n", famp->seq, famp->fa->id, famp->fa->seq, famp->mo->id, famp->mo->seq);*/
        if (famp->fa->fam) {
            addLink(linkList, nlink, linkInd, famp->seq, famp->fa->fam->seq,
                    famp->fa->seq);
            addLink(linkList, nlink, linkInd, famp->fa->fam->seq, famp->seq,
                    famp->fa->seq);
        }
        if (famp->mo->fam) {
            addLink(linkList, nlink, linkInd, famp->seq, famp->mo->fam->seq,
                    famp->mo->seq);
            addLink(linkList, nlink, linkInd, famp->mo->fam->seq, famp->seq,
                    famp->mo->seq);
        }
        famp2 = pedp->fam1;
        while (famp2 != famp) {
            if (famp2->fa == famp->fa) {
                addLink(linkList, nlink, linkInd, famp->seq, famp2->seq,
                        famp2->fa->seq);
                addLink(linkList, nlink, linkInd, famp2->seq, famp->seq,
                        famp2->fa->seq);
            }
            if (famp2->mo == famp->mo) {
                addLink(linkList, nlink, linkInd, famp->seq, famp2->seq,
                        famp2->mo->seq);
                addLink(linkList, nlink, linkInd, famp2->seq, famp->seq,
                        famp2->mo->seq);
            }
            famp2 = famp2->next;
        }
        famp = famp->next;
    }

    do {
        done = TRUE;
        for (i = 0; i < NumFam; i++) {
            if (nlink[i] == 1) {
                for (j = 0; j < NumFam; j++)
                    rmLink(linkList, nlink, linkInd, (int) j, (int) i);
                linkp = linkList[i];
                while (linkp) {
                    linkInd[linkp->ind]--;
                    lastp = linkp;
                    linkp = linkp->next;
                    free(lastp);
                }
                nlink[i] = 0;
                linkList[i] = NULL;
                done = FALSE;
            }
/*else {
linkp = linkList[i];
while (linkp) {
printf("fam1 = %d  fam2 = %d  id = %d\n", i, linkp->fam, linkp->ind);
linkp = linkp->next;
}
}*/
        }
    } while (!done);
}

int findBreaks(int *nlink, int *linkInd, char *lbrkId)
{
    int i, nlbrk, nodes, narcs;
    Ind *indp;

    narcs = 0;
    nodes = 0;

    for (i = 0; i < NumFam; i++) {
        if (nlink[i]) {
            narcs += nlink[i];
            nodes++;
        }
    }

    for (i = 0; i < NumInd; i++)
        if (linkInd[i])
            nodes++;

    if (narcs >= nodes)
        nlbrk = narcs - nodes + 1;

    if (nlbrk > 1)
        return nlbrk;

    *lbrkId = '\0';
    for (i = 0; !*lbrkId && i < NumInd; i++) {
        if (linkInd[i]) {
            indp = IndArray[IndSort[i]];
            if (indp->fam && nlink[indp->fam->seq])
                strcpy(lbrkId, indp->id);
        }
    }

    return nlbrk;
}

void addLink(Link **linkList, int *nlink, int *linkInd, int fam1, int fam2,
             int ind)
{
    int found;
    Link *linkp, *lastp;

    found = FALSE;
    linkp = linkList[fam1];
    lastp = NULL;
    while (linkp) {
        if (linkp->ind == ind)
            found = TRUE;
        lastp = linkp;
        linkp = linkp->next;
    }

    linkp = (Link *) allocMem(sizeof(Link));
    linkp->ind = ind;
    linkp->fam = fam2;
    linkp->next = NULL;
    linkInd[ind]++;

    if (lastp)
        lastp->next = linkp;
    else
        linkList[fam1] = linkp;

    if (!found)
        nlink[fam1]++;
}

void rmLink(Link **linkList, int *nlink, int *linkInd, int fam1, int fam2)
{
    int ind;
    Link *linkp, *lastp, *nextp;

    if (!nlink[fam1])
        return;

    linkp = linkList[fam1];
    lastp = NULL;
    while (linkp) {
        nextp = linkp->next;
        if (linkp->fam == fam2) {
            ind = linkp->ind;
            linkInd[ind]--;
            if (lastp)
                lastp->next = nextp;
            else
                linkList[fam1] = nextp;
            free(linkp);
            linkp = linkList[fam1];
            while (linkp) {
                if (linkp->ind == ind)
                    return;
                linkp = linkp->next;
            }
            nlink[fam1]--;
            if (!nlink[fam1])
                linkList[fam1] = NULL;
            return;
        }
        lastp = linkp;
        linkp = linkp->next;
    }
}

void warshall (unsigned char **adj, int n)
{
    int i, j, k;
    unsigned char *adji, *adjj;

    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            if (adj[i][j])
                for (k = 0, adji = adj[i], adjj = adj[j];
                        k < n; k++, adji++, adjj++)
                    *adji = *adji || *adjj;
}

void assignSeq (void)
{
    int i, w1, w2, width;
    int famseq, curped;
    int indseq; 
    char fmt[1024], *seqList[MAXIND];

    w1 = log((double)NumPed) / log(10.) + 1;
    w2 = log((double)NumInd) / log(10.) + 1;
    sprintf(fmt, "%%%dd%%%dd%%%dd%%%dd", w1, w2, w2, w2);
    width = w1 + 3*w2;
    for (i = 0; i < NumInd; i++) {
        seqList[i] = (char *) allocMem((size_t)(width+1));
        if (IndArray[i]->fam)
            famseq = IndArray[i]->fam->seq;
        else
            famseq = 0;
        indseq = IndArray[findInd(IndArray[i]->id)]->seq;
        sprintf(seqList[i], fmt, IndArray[i]->ped, IndArray[i]->gen,
                famseq, indseq);
    }

    qSort(seqList, width, (int) NumInd, IndSeq, FALSE);

    curped = IndArray[IndSeq[0]]->ped;
    PedArray[curped]->seq1 = 0;
    for (i = 0; i < NumInd; i++) {
        if (IndArray[IndSeq[i]]->ped != curped) {
            curped = IndArray[IndSeq[i]]->ped;
            PedArray[curped]->seq1 = (int) i;
        }
        IndArray[IndSeq[i]]->seq = (int) i;
        free(seqList[i]);
    }
}

void calcKin2 (void)
{
    int i, j, ifa, imo, jfa, jmo, count;
    int n, itwin[MAXIND], twin1[MXTWIN];
    int itwinid;
    int rc, system();
    float *kin2[MAXIND], delta7;
    FILE *outfp;

    for (i = 0; i < MXTWIN; i++)
        twin1[i] = -1;

    n = 0;
    for (i = 0; i < NumInd; i++) {
        itwin[i] = i;
        itwinid = IndArray[IndSeq[i]]->itwinid;
        if (itwinid) {
            if (twin1[itwinid-1] != -1)
                itwin[i] = twin1[itwinid-1];
            else
                twin1[itwinid-1] = i;
        }
    }

    count = 0;
    for (i = 0; i < NumInd; i++) {
        if (itwin[i] == i) n++;
        kin2[i] = (float *) allocMem((size_t)((i+1)*sizeof(float)));
        for (j = 0; j < i; j++)
            kin2[i][j] = 0;
        if (IndArray[IndSeq[i]]->fam == NULL) {
            count++;
            kin2[i][i] = 1;
        }
        else
            kin2[i][i] = 0;
    }

    do {
        for (i = 0; i < NumInd; i++) {
            if (itwin[i] != i || kin2[i][i] != 0) continue;
            if (IndArray[IndSeq[i]]->fam == NULL) continue;
            ifa = itwin[IndArray[IndSeq[i]]->fam->fa->seq];
            imo = itwin[IndArray[IndSeq[i]]->fam->mo->seq];
            if (kin2[ifa][ifa] == 0 || kin2[imo][imo] == 0) continue;
            for (j = 0; j < NumInd; j++) {
                if (itwin[j] != j || kin2[j][j] == 0) continue;
                kin2[max(i,j)][min(i,j)] =
                                 .5 * ( kin2[max(ifa,j)][min(ifa,j)] +
                                        kin2[max(imo,j)][min(imo,j)] );
            }
            count++;
            kin2[i][i] = 1 + .5 * kin2[max(ifa,imo)][min(ifa,imo)];
        }
    } while (count < n);

    for (i = 0; i < NumInd; i++) {
        for (j = 0; j < i; j++)
            kin2[i][j] = kin2[max(itwin[i],itwin[j])][min(itwin[i],itwin[j])];
        kin2[i][i] = kin2[itwin[i]][itwin[i]];
    }

    isInbred = FALSE;
    for (i = 0; i < NumPed; i++)
        PedArray[i]->inbred = FALSE;

    outfp = openFile("phi2", "w");

    for (i = 0; i < NumInd; i++) {
        for (j = 0; j < i; j++) {
            if (IndArray[IndSeq[i]]->ped != IndArray[IndSeq[j]]->ped)
                continue;

            if (itwin[i] == itwin[j])
                delta7 = 1;

            else {
                delta7 = 0;
                if (IndArray[IndSeq[i]]->fam != NULL &&
                        IndArray[IndSeq[j]]->fam != NULL) {
                    ifa = IndArray[IndSeq[i]]->fam->fa->seq;
                    imo = IndArray[IndSeq[i]]->fam->mo->seq;
                    jfa = IndArray[IndSeq[j]]->fam->fa->seq;
                    jmo = IndArray[IndSeq[j]]->fam->mo->seq;
                    delta7 = .25 * ( kin2[max(ifa,jfa)][min(ifa,jfa)] *
                                     kin2[max(imo,jmo)][min(imo,jmo)] +
                                     kin2[max(ifa,jmo)][min(ifa,jmo)] *
                                     kin2[max(imo,jfa)][min(imo,jfa)] );
                }
            }

            if (kin2[i][j])
                fprintf(outfp, "%8d %8d %10.7f %10.7f\n", i + 1, j + 1,
                        kin2[i][j], delta7);
        }

        fprintf(outfp, "%8d %8d %10.7f %10.7f\n", i + 1, i + 1,
                kin2[i][i], 1.);

        if (kin2[i][i] > 1.) {
            isInbred = TRUE;
            PedArray[IndArray[IndSeq[i]]->ped]->inbred = TRUE;
        }
    }

    fclose(outfp);
    rc = system("gzip -f phi2");

    for (i = 0; i < NumInd; i++)
        free(kin2[i]);
}

void makeHHoldMat (void)
{
    int i, j;
    int rc, system();
    FILE *outfp;

    outfp = openFile("house", "w");

    for (i = 0; i < NumInd; i++) {
        for (j = 0; j < i; j++) {
            if (strlen(IndArray[IndSeq[i]]->hhid) &&
                    !strcmp(IndArray[IndSeq[i]]->hhid,
                            IndArray[IndSeq[j]]->hhid))
                fprintf(outfp, "%5d %5d %10.7f %10.7f\n", i + 1, j + 1,
                        1., 0.);
        }
        fprintf(outfp, "%5d %5d %10.7f %10.7f\n", i + 1, i + 1, 1., 0.);
    }

    fclose(outfp);
    rc = system("gzip -f house");
}

void writeIndex (void)
{
    int i, done;
    int iseq;
    FILE *outfp;

    outfp = openFile("pedindex.out", "w");

    iseq = 0;
    for (i = 0; i < NumPed; i++) {
        done = FALSE;
        while (iseq < NumInd && !done) {
            if (IndArray[IndSeq[iseq]]->ped == (int) i) {
                if (IndArray[IndSeq[iseq]]->fam)
		    if (TwinOutLen <= 3)
		    {
                    fprintf(outfp, "%8d %8d %8d %1d %3d %8d %8d %s\n",
                            IndArray[IndSeq[iseq]]->seq + 1,
                            IndArray[IndSeq[iseq]]->fam->fa->seq + 1,
                            IndArray[IndSeq[iseq]]->fam->mo->seq + 1,
                            IndArray[IndSeq[iseq]]->sex,
                            IndArray[IndSeq[iseq]]->itwinid,
                            IndArray[IndSeq[iseq]]->ped + 1,
                            IndArray[IndSeq[iseq]]->gen,
                            IndArray[IndSeq[iseq]]->id);
		    } else {
                    fprintf(outfp, "%8d %8d %8d %1d %5d %8d %8d %s\n",
                            IndArray[IndSeq[iseq]]->seq + 1,
                            IndArray[IndSeq[iseq]]->fam->fa->seq + 1,
                            IndArray[IndSeq[iseq]]->fam->mo->seq + 1,
                            IndArray[IndSeq[iseq]]->sex,
                            IndArray[IndSeq[iseq]]->itwinid,
                            IndArray[IndSeq[iseq]]->ped + 1,
                            IndArray[IndSeq[iseq]]->gen,
                            IndArray[IndSeq[iseq]]->id);
		    }
                else
		    if (TwinOutLen <= 3)
		    {
                    fprintf(outfp, "%8d %8d %8d %1d %3d %8d %8d %s\n",
                            IndArray[IndSeq[iseq]]->seq + 1,
                            0, 0,
                            IndArray[IndSeq[iseq]]->sex,
                            IndArray[IndSeq[iseq]]->itwinid,
                            IndArray[IndSeq[iseq]]->ped + 1,
                            IndArray[IndSeq[iseq]]->gen,
                            IndArray[IndSeq[iseq]]->id);
		    } else {
                    fprintf(outfp, "%8d %8d %8d %1d %8d %8d %8d %s\n",
                            IndArray[IndSeq[iseq]]->seq + 1,
                            0, 0,
                            IndArray[IndSeq[iseq]]->sex,
                            IndArray[IndSeq[iseq]]->itwinid,
                            IndArray[IndSeq[iseq]]->ped + 1,
                            IndArray[IndSeq[iseq]]->gen,
                            IndArray[IndSeq[iseq]]->id);
		    }
                iseq++;
            }
            else if (IndArray[IndSeq[iseq]]->ped == -1)
                iseq++;
            else
                done = TRUE;
        }
    }

    fclose(outfp);

    outfp = openFile("pedindex.cde", "w");

    fprintf(outfp,
            "pedindex.out                                          \n");
    fprintf(outfp,
            " 8 IBDID                 IBDID                       I\n");
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");
    fprintf(outfp,
            " 8 FATHER'S IBDID        FIBDID                      I\n");
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");
    fprintf(outfp,
            " 8 MOTHER'S IBDID        MIBDID                      I\n");
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");
    fprintf(outfp,
            " 1 SEX                   SEX                         I\n");
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");
    if (TwinOutLen > 3)
    {
	fprintf(outfp,
            " 8 MZTWIN                MZTWIN                      I\n");
    } else {
	fprintf(outfp,
            " 3 MZTWIN                MZTWIN                      I\n");
    }
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");
    fprintf(outfp,
            " 8 PEDIGREE NUMBER       PEDNO                       I\n");
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");
    fprintf(outfp,
            " 8 GENERATION NUMBER     GEN                         I\n");
    fprintf(outfp,
            " 1 BLANK                 BLANK                       C\n");

    if (FamidLen) {
        fprintf(outfp,
                "%2d FAMILY ID             FAMID                       C\n",
                FamidLen);
    }
    fprintf(outfp,
            "%2d ID                    ID                          C\n",
            IdLen);

    fclose(outfp);
}

void writeMCarloFiles (int loc)
{
    int i, j, done;
    int iseq, all1, all2;
    Ind *indp;
    Loc *locp = LocArray[loc];
    char dirname[1024], outfile[1024];
    char twinid[4];
    FILE *outfp;

    sprintf(dirname, "d_%s", locp->mrkName);
    makeDir(dirname, 0755);

    sprintf(outfile, "%s/translat.tab", dirname);
    outfp = openFile(outfile, "w");

    fprintf(outfp,"(I6,2X,A8)\n");
    fprintf(outfp,"(3A5,A1,A3,A6)\n");

    iseq = 0;
    for (i = 0; i < NumPed; i++) {
        fprintf(outfp, "%6d  FAM%05d\n", PedArray[i]->nind, i + 1);
        done = FALSE;
        while (iseq < NumInd && !done) {
            indp = IndArray[iseq];
            if (indp->itwinid)
               sprintf(twinid, "%3d", indp->itwinid);
            else
               sprintf(twinid, "   ");
            if (indp->ped == (int) i) {
                if (indp->fam) {
                    if (indp->mrkall[loc][0] != -1) {
                        all1 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][0]]) + 1;
                        all2 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][1]]) + 1;
                        fprintf(outfp, "%5d%5d%5d%d%s%3d%3d\n",
                                indp->seq + 1, indp->fam->fa->seq + 1,
                                indp->fam->mo->seq + 1, indp->sex, twinid,
                                min(all1,all2), max(all1,all2));
                    }
                    else
                        fprintf(outfp, "%5d%5d%5d%d%s      \n",
                                indp->seq + 1, indp->fam->fa->seq + 1,
                                indp->fam->mo->seq + 1, indp->sex, twinid);
                }
                else {
                    if (indp->mrkall[loc][0] != -1) {
                        all1 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][0]]) + 1;
                        all2 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][1]]) + 1;
                        fprintf(outfp, "%5d          %d%s%3d%3d\n",
                                indp->seq + 1, indp->sex, twinid,
                                min(all1,all2), max(all1,all2));
                    }
                    else
                        fprintf(outfp, "%5d          %d%s      \n",
                                indp->seq + 1, indp->sex, twinid);
                }
                iseq++;
            }
            else if (indp->ped == -1)
                iseq++;
            else
                done = TRUE;
        }
    }

    fclose(outfp);

    sprintf(outfile, "%s/ibd.loc", dirname);
    outfp = openFile(outfile, "w");

    if (xLinked)
        fprintf(outfp,"%-8.8sX-LINKED%2d%3d\n", locp->mrkName, locp->numAll,
                locp->numAll*(locp->numAll+1)/2);
    else
        fprintf(outfp,"%-8.8sAUTOSOME%2d%3d\n", locp->mrkName, locp->numAll,
                locp->numAll*(locp->numAll+1)/2);

    for (i = 0; i < locp->numAll; i++)
        fprintf(outfp, "%2d      %8.7f\n", i + 1,
                       locp->allFreq[locp->allSort[i]]);

    for (i = 0; i < locp->numAll; i++) {
        for (j = i; j < locp->numAll; j++) {
            fprintf(outfp, " %2d %2d   1\n", i + 1, j + 1);
            fprintf(outfp, "%2d/%2d\n", i + 1, j + 1);
        }
    }

    fclose(outfp);

    sprintf(outfile, "%s/ibd.bat", dirname);
    outfp = openFile(outfile, "w");

    fprintf(outfp,"9\n");
    fprintf(outfp,"%-8.8s\n", locp->mrkName);
    fprintf(outfp,"\n");
    fprintf(outfp,"21\n");
    fprintf(outfp,"n\n");
    fclose(outfp);
}

void writeMLEfreqFiles (int loc)
{
    int i, done;
    int untyped[MAXPED];
    int iseq, all1, all2;
    Ind *indp;
    Loc *locp = LocArray[loc];
    char dirname[1024], outfile[1024];
    char twinid[4];
    FILE *outfp;

    sprintf(dirname, "d_%s", locp->mrkName);
    makeDir(dirname, 0755);

    sprintf(outfile, "%s/allfreq.ped", dirname);
    outfp = openFile(outfile, "w");

    fprintf(outfp,"(I6,2X,A8)\n");
    fprintf(outfp,"(3A5,A1,A3,A5)\n");

    iseq = 0;
    for (i = 0; i < NumPed; i++) {
        untyped[i] = TRUE;
        done = FALSE;
        while (iseq < NumInd && !done) {
            indp = IndArray[iseq];
            if (indp->ped == (int) i) {
                if (indp->mrkall[loc][0] != -1) untyped[i] = FALSE;
                iseq++;
            }
            else if (indp->ped == -1)
                iseq++;
            else
                done = TRUE;
        }
    }

    iseq = 0;
    for (i = 0; i < NumPed; i++) {

        /* skip pedigrees in which no individual is typed */
        if (untyped[i]) {
            done = FALSE;
            while (iseq < NumInd && !done) {
                indp = IndArray[iseq];
                if (indp->ped == (int) i || indp->ped == -1)
                    iseq++;
                else
                    done = TRUE;
            }
            continue;
        }

        fprintf(outfp, "%6d  FAM%05d\n", PedArray[i]->nind, i + 1);
        done = FALSE;
        while (iseq < NumInd && !done) {
            indp = IndArray[iseq];
            if (indp->itwinid)
                sprintf(twinid, "%3d", indp->itwinid);
            else
                sprintf(twinid, "   ");
            if (indp->ped == (int) i) {
                if (indp->fam) {
                    if (indp->mrkall[loc][0] != -1) {
                        all1 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][0]]) + 1;
                        all2 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][1]]) + 1;
                        fprintf(outfp, "%5d%5d%5d%d%s%2d/%2d\n",
                                indp->seq + 1, indp->fam->fa->seq + 1,
                                indp->fam->mo->seq + 1, indp->sex, twinid,
                                all1, all2);
                    }
                    else
                        fprintf(outfp, "%5d%5d%5d%d%s     \n",
                                indp->seq + 1, indp->fam->fa->seq + 1,
                                indp->fam->mo->seq + 1, indp->sex, twinid);
                }
                else {
                    if (indp->mrkall[loc][0] != -1) {
                        all1 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][0]]) + 1;
                        all2 = findAllele((int)loc,
                                   locp->allList[indp->mrkall[loc][1]]) + 1;
                        fprintf(outfp, "%5d          %d%s%2d/%2d\n",
                                indp->seq + 1, indp->sex, twinid, all1, all2);
                    }
                    else
                        fprintf(outfp, "%5d          %d%s     \n",
                                indp->seq + 1, indp->sex, twinid);
                }
                iseq++;
            }
            else if (indp->ped == -1)
                iseq++;
            else
                done = TRUE;
        }
    }

    fclose(outfp);

    sprintf(outfile, "%s/allfreq.loc", dirname);
    outfp = openFile(outfile, "w");

    if (xLinked)
        fprintf(outfp,"%-8.8sX-LINKED%2d\n", locp->mrkName, locp->numAll);
    else
        fprintf(outfp,"%-8.8sAUTOSOME%2d\n", locp->mrkName, locp->numAll);

    for (i = 0; i < locp->numAll; i++)
        fprintf(outfp, "%5d   %8.7f\n", i + 1,
                locp->allFreq[locp->allSort[i]]);

    fclose(outfp);

    sprintf(outfile, "%s/allfreq.bat", dirname);
    outfp = openFile(outfile, "w");

    fprintf(outfp,"9\n");
    fprintf(outfp,"%-8.8s\n", locp->mrkName);
    fprintf(outfp,"17\n");
    fprintf(outfp,"%2d\n", locp->numAll);
    fprintf(outfp,"21\n");
    fprintf(outfp,"n\n");
    fclose(outfp);

    sprintf(outfile, "%s/allfreq.mod", dirname);
    outfp = openFile(outfile, "w");

    for (i = 0; i < locp->numAll; i++)
        fprintf(outfp,
                "%2d%-5s      %8.6fD+00   0.100000D-05   0.100000D+01\n",
                i + 1, locp->allList[locp->allSort[i]],
                locp->allFreq[locp->allSort[i]]);
    fprintf(outfp, "CNS LINES=%2d\n", locp->numAll);
    for (i = 0; i < locp->numAll; i++)
        fprintf(outfp, "  1 %2d 0.1D+01\n", i + 1);
    fprintf(outfp, "CVALUES  = 1\n");
    fprintf(outfp, "     1 0.1D+01\n");
    fclose(outfp);
}

void writeLinkageFiles (int loc)
{
    int i;
    int all1, all2;
    Ind *indp;
    Loc *locp = LocArray[loc];
    char dirname[1024], outfile[1024];
    FILE *outfp;

    sprintf(dirname, "d_%s", locp->mrkName);
    makeDir(dirname, 0755);

    sprintf(outfile, "%s/ped.raw", dirname);
    outfp = openFile(outfile, "w");

    for (i = 0; i < NumInd; i++) {
        indp = IndArray[i];
        if (dropSingles && PedArray[indp->ped]->nind == 1)
            continue;
        if (indp->mrkall[loc][0] == -1)
            all1 = all2 = 0;
        else {
            all1 = findAllele((int)loc,
                              locp->allList[indp->mrkall[loc][0]]) + 1;
            all2 = findAllele((int)loc,
                              locp->allList[indp->mrkall[loc][1]]) + 1;
        }
        if (indp->fam)
            fprintf(outfp, "%5d %5d %5d %5d %1d %1d %3d %2d %2d\n",
                    indp->ped + 1, indp->seq + 1, indp->fam->fa->seq + 1,
                    indp->fam->mo->seq + 1, indp->sex, 0, indp->itwinid,
                    all1, all2);
        else
            fprintf(outfp, "%5d %5d %5d %5d %1d %1d %3d %2d %2d\n",
                    indp->ped + 1, indp->seq + 1, 0, 0, indp->sex, 0,
                    indp->itwinid, all1, all2);
    }

    fclose(outfp);

    sprintf(outfile, "%s/datafile.dat", dirname);
    outfp = openFile(outfile, "w");

    if (xLinked)
        fprintf(outfp, "2 1 1 5\n");
    else
        fprintf(outfp, "2 1 0 5\n");
    fprintf(outfp, "0 0.00000000 0.00000000 0\n");
    fprintf(outfp, " 1 2\n\n");
    fprintf(outfp, "1 2\n");
    fprintf(outfp, " 0.99999999 0.00000001\n");
    fprintf(outfp, " 1\n");
    fprintf(outfp, " 0.00000000 0.00000000 1.00000000\n");
    if (xLinked)
        fprintf(outfp, " 0.00000000 0.50000000\n");
    fprintf(outfp, "2\n\n");

    if (locp->numAll > 1) {
        fprintf(outfp,"3 %d\n", locp->numAll);
        fprintf(outfp,"%11.8f", locp->allFreq[locp->allSort[0]]);
        for (i = 1; i < locp->numAll; i++)
            fprintf(outfp, "%11.8f", locp->allFreq[locp->allSort[i]]);
        fprintf(outfp, "\n");
    }
    else {
        fprintf(outfp,"3 2\n");
        fprintf(outfp," 0.90000000 0.10000000\n");
    }

    fprintf(outfp, "\n0 0\n");
    fprintf(outfp, " 0.00000000\n");
    fprintf(outfp, "1 0.10000000 0.09000000\n");

    fclose(outfp);
}

void writeMakepedCmd (int loc)
{
    int i;
    Loc *locp = LocArray[loc];
    char dirname[1024], outfile[1024];
    FILE *outfp;

    sprintf(dirname, "d_%s", locp->mrkName);
    sprintf(outfile, "%s/makeped.cmd", dirname);
    outfp = openFile(outfile, "w");

    fprintf(outfp, "ped.raw\n");
    fprintf(outfp, "pedin.dat\n");

    if (MaxLbrk) {
        fprintf(outfp, "y\n");
        fprintf(outfp, "n\n");
        for (i = 0; i < NumPed; i++) {
            if (PedArray[i]->hasloops)
                fprintf(outfp, "%d\n%d\n", i + 1,
                        IndArray[PedArray[i]->lbrkind]->seq + 1);
        }
        fprintf(outfp, "0\n");
        fprintf(outfp, "n\n");
        fprintf(outfp, "y\n");
    }
    else {
        fprintf(outfp, "n\n");
        fprintf(outfp, "y\n");
    }

    fclose(outfp);
}

void writeMMSibsFiles (void)
{
    int i, loc;
    int all1, all2;
    char chrnum[1024];
    float mrkloc[MAXLOC];
    char mrknam[1024];
    Ind *indp;
    Loc *locp;
    FILE *infp, *outfp;

    outfp = openFile("sibs.ped", "w");

    for (i = 0; i < NumInd; i++) {
        indp = IndArray[i];
        if (dropSingles && PedArray[indp->ped]->nind == 1)
            continue;

        if (indp->fam)
            fprintf(outfp, "%6d%6s%6s%6s%3d%3d", indp->ped + 1, indp->id,
                    indp->fam->fa->id, indp->fam->mo->id, indp->sex, 2);
        else
            fprintf(outfp, "%6d%6s%6d%6d%3d%3d", indp->ped + 1, indp->id,
                    0, 0, indp->sex, 1);

        locp = LocArray[0];
        if (indp->mrkall[0][0] == -1)
            all1 = 0;
        else
            all1 = findAllele((int)0,
                              locp->allList[indp->mrkall[0][0]]) + 1;
        if (indp->mrkall[0][1] == -1)
            all2 = 0;
        else {
            all2 = findAllele((int)0,
                              locp->allList[indp->mrkall[0][1]]) + 1;
            if (!all1) all1 = all2;
        }
        fprintf(outfp, "%3d%3d", all1, all2);
        for (loc = 1; loc < NumLoc; loc++) {
            locp = LocArray[loc];
            if (indp->mrkall[loc][0] == -1)
                all1 = 0;
            else
                all1 = findAllele((int)loc,
                                  locp->allList[indp->mrkall[loc][0]]) + 1;
            if (indp->mrkall[loc][1] == -1)
                all2 = 0;
            else {
                all2 = findAllele((int)loc,
                                  locp->allList[indp->mrkall[loc][1]]) + 1;
                if (!all1) all1 = all2;
            }
            fprintf(outfp, "  %3d%3d", all1, all2);
        }
        fprintf(outfp, "\n");
    }

    fclose(outfp);

    infp = openFile(MapFile, "r");
    if (fscanf(infp, "%s\n", chrnum) != 1) {
        sprintf(ErrMsg, "invalid record, line 1 of map-data file");
        fatalError();
    }
    for (loc = 0; loc < NumLoc; loc++) {
        if (fscanf(infp, "%s %f\n", mrknam, &mrkloc[loc]) != 2) {
            sprintf(ErrMsg, "invalid record, line %d of map-data file",
                    loc + 2);
            fatalError();
        }
    }
    fclose(infp);

    outfp = openFile("sibs.loc", "w");

    if (xLinked)
        fprintf(outfp, "%2d 1 1 5\n", NumLoc + 1);
    else
        fprintf(outfp, "%2d 1 0 5\n", NumLoc + 1);
    fprintf(outfp, "0 0.0 0.0 0\n");
    fprintf(outfp, "1");

    for (i = 2; i < NumLoc + 2; i++)
        fprintf(outfp, "%3d", i);
    fprintf(outfp, "\n");

    fprintf(outfp, "1 2\n");
    fprintf(outfp, "0.990000 0.100000\n");
    fprintf(outfp, "1\n");
    fprintf(outfp, "0.001000 0.001000 0.999000\n");
    if (xLinked)
        fprintf(outfp, "0.001000 0.499000\n");

    for (loc = 0; loc < NumLoc; loc++) {
        locp = LocArray[loc];
        fprintf(outfp, "3        %2d\n", locp->numAll);
        fprintf(outfp, "%8.6f", locp->allFreq[locp->allSort[0]]);
        for (i = 1; i < locp->numAll; i++)
            fprintf(outfp, " %8.6f", locp->allFreq[locp->allSort[i]]);
        fprintf(outfp, "\n");
    }

    fprintf(outfp, "0 0\n");
    fprintf(outfp, "%5.1f", mrkloc[0]);
    for (loc = 1; loc < NumLoc; loc++)
        fprintf(outfp, " %5.1f", mrkloc[loc] - mrkloc[loc-1]);
    fprintf(outfp, "\n");
    fprintf(outfp, "1 0 0.5\n");

    fclose(outfp);
}

/*
 *  allSort[ findAllele( allList[ n ] ) ] = n
 */
int findAllele (int loc, char *allele)
{
    int ndx, lo, hi, cmp, a0, a1;
    Loc *locp = LocArray[loc];

    lo = 0;
    hi = locp->numAll - 1;
    while (hi >= lo) {
        ndx = (lo + hi)/2;
        if (locp->allNumeric) {
            a0 = atoi(allele);
            a1 = atoi(locp->allList[locp->allSort[ndx]]);
            if (a0 < a1)
                cmp = -1;
            else if (a0 > a1)
                cmp = 1;
            else
                cmp = 0;
        }
        else
            cmp = strcmp(allele, locp->allList[locp->allSort[ndx]]);
        if (cmp == 0)
            return ndx;
        else if (cmp < 0)
            hi = ndx - 1;
        else
            lo = ndx + 1;
    }

    return -1;
}

void writeInfo (void)
{
    int i;
    FILE *outfp;

    if (doIndex) {
        outfp = fopen("pedigree.info", "a");
        if (!outfp) {
            sprintf(ErrMsg, "cannot open pedigree.info");
            fatalError();
        }
        fprintf(outfp, "%d %d %d %d %d\n", IdLen, SexLen, TwinIdLen, HHIdLen,
                FamidLen);

        /* add singletons to count of nuclear families */
        for (i = 0; i < NumPed; i++) {
            if (PedArray[i]->nfou == 1) {
                PedArray[i]->nfam = 1;
                NumFam++;
            }
        }
        fprintf(outfp, "%d %d %d %d\n", NumPed, NumFam, NumInd, NumFou);
        for (i = 0; i < NumPed; i++)
            fprintf(outfp, "%d %d %d %d %c\n", PedArray[i]->nfam,
                    PedArray[i]->nind, PedArray[i]->nfou,
                    PedArray[i]->nlbrk, PedArray[i]->inbred ? 'y' : 'n');
        fclose(outfp);
    }
    else {
        outfp = fopen("marker.info", "a");
        if (!outfp) {
            sprintf(ErrMsg, "cannot open marker.info");
            fatalError();
        }
        for (i = 0; i < NumLoc; i++)
            fprintf(outfp, "%s %d %d\n", LocArray[i]->mrkName,
                    LocArray[i]->numTyp, LocArray[i]->numFouTyp);
        fclose(outfp);
    }
}

int getAlleles (char *gtype, char **allele, int *allnum)
{
    int i;
    int divided = FALSE, numeric = FALSE;
    char *p, *q, *start;

    p = gtype;
    *allele[0] = '\0';
    *allele[1] = '\0';

    while (*p) {
        if (*p == '(' || *p == ')') *p = ' ';
        p++;
    }

    p = gtype;
    while (*p && (*p == ' ' || *p == '\t')) p++;
    if (!*p) {
        return TRUE;
    }

    start = p;
    while (*p && *p != '/' && *p != ' ' && *p != '\t') p++;
    while (*p && (*p == ' ' || *p == '\t')) p++;
    if (*p)
        divided = TRUE;

    p = start;
    if (!divided) {
        if (isdigit(*p))
            numeric = TRUE;
        else if (!isalpha(*p))
            return FALSE;
    }

    for (i = 0; i < 2; i++) {
        q = allele[i];
        *q = *p;

        if (divided) {
            while (*p && *p != '/' && *p != ' ' && *p != '\t')
                *q++ = *p++;
            *q = '\0';
            if (!strlen(allele[i]) && !xLinked)
                return FALSE;
            while (*p && (*p == ' ' || *p == '\t' || *p == '/')) p++;
        }
        else if (numeric) {
            if (!isdigit(*p))
                return FALSE;
            *q++ = *p++;
            while (*p && (isalpha(*p) || *p == '\''))
                *q++ = *p++;
            *q = '\0';
            if (!strlen(allele[i]))
                return FALSE;
        }
        else {
            if (!isalpha(*p))
                return FALSE;
            *q++ = *p++;
            while (*p && (isdigit(*p) || *p == '\''))
                *q++ = *p++;
            *q = '\0';
            if (!strlen(allele[i]))
                return FALSE;
        }
    }

    while (*p && (*p == ' ' || *p == '\t')) p++;
    if (*p)
        return FALSE;

    for (i = 0; i < 2; i++) {
        q = allele[i];
        while (*q) {
            *allnum = *allnum && isdigit(*q);
            q++;
        }
    }

    if (*allele[0] && !strcmp(allele[0], "0")) *allele[0] = '\0';
    if (*allele[1] && !strcmp(allele[1], "0")) *allele[1] = '\0';

    if (*allele[0] && !strcmp(allele[0], "-")) *allele[0] = '\0';
    if (*allele[1] && !strcmp(allele[1], "-")) *allele[1] = '\0';

    if (!xLinked &&
            (!*allele[0] && *allele[1] || *allele[0] && !*allele[1]))
        return FALSE;

    return TRUE;
}

void cntAlleles (int loc, char **allele, int *allCnt, int *mrkall)
{
    int i, j, found;
    int temp;
    Loc *locp = LocArray[loc];

    if (!strlen(allele[1])) {
        mrkall[0] = -1;
        mrkall[1] = -1;
        return;
    }

    i = 0;
    if (!strlen(allele[0])) {
        mrkall[0] = -1;
        i = 1;
    }
    for ( ; i < 2; i++) {
        found = FALSE;
        for (j = 0; !found && j < locp->numAll; j++) {
            if (!strcmp(allele[i], locp->allList[j])) {
                found = TRUE;
                mrkall[i] = j;
                if (locp->noLocInfo) allCnt[j]++;
            }
        }
        if (!found) {
            if (!locp->noLocInfo) {
                sprintf(ErrMsg, "unknown allele [%s] found for marker %s",
                        allele[i], locp->mrkName);
                fatalError();
            }
            if (locp->numAll == MAXALL) {
                sprintf(ErrMsg, "locus %d has too many alleles, MAXALL = %d",
                        loc + 1, MAXALL);
                fatalError();
            }
            locp->allList[locp->numAll] =
                          (char *) allocMem((size_t)(GtypeLen+1));
            strcpy(locp->allList[locp->numAll], allele[i]);
            allCnt[locp->numAll] = 1;
            mrkall[i] = locp->numAll;
            locp->numAll++;
        }
    }

    if (mrkall[0] >= mrkall[1]) {
        temp = mrkall[0];
        mrkall[0] = mrkall[1];
        mrkall[1] = temp;
    }
}

void qSort (char **vals, int vlen, int nvals, int *ord, int numeric)
{
    int i, j, t, l, h;
    int ip;
    char *mid;
 
#define STKSIZE 1000
    int lo[STKSIZE], hi[STKSIZE];

    if (!nvals || !vlen)
        return;

    mid = (char *) allocMem((size_t)(vlen+1));

    for (i = 0; i < nvals; i++)
        ord[i] = i;

    if (nvals == 1)
        return;

    t = 0;
    lo[0] = 0;
    hi[0] = nvals - 1;
    do {
        l = lo[t];
        h = hi[t];
        t--;
        do {
            i = l;
            j = h;
            strcpy(mid, vals[ord[(i+j)/2+1]]);
            do {
                if (numeric) {
                    while (atoi(vals[ord[i]]) < atoi(mid)) i++;
                    while (atoi(mid) < atoi(vals[ord[j]])) j--;
                }
                else {
                    while (strcmp(vals[ord[i]], mid) < 0) i++;
                    while (strcmp(mid, vals[ord[j]]) < 0) j--;
                }
                if (i <= j) {
                    ip = ord[i];
                    ord[i] = ord[j];
                    ord[j] = ip;
                    i++;
                    j--;
                }
            } while (i <= j);
            if (i < h) {
                t++;
                if (t >= STKSIZE) {
                    strcpy(ErrMsg, "qSort: stack size exceeded");
                    fatalError();
                }
                lo[t] = i;
                hi[t] = h;
            }
            h = j;
        } while (l < h);
    } while (t >= 0);
}

void *allocMem (size_t nbytes)
{
    void *ptr;
    ptr = (void *) malloc(nbytes);
    if (!ptr) {
        strcpy(ErrMsg, "not enough memory");
        fatalError();
    }
    return ptr;
}

FILE *openFile (char *fileName, char *fileMode)
{
    FILE *fp;
    fp = fopen(fileName, fileMode);
    if (!fp) {
        sprintf(ErrMsg, "cannot open file \"%s\"", fileName);
        fatalError();
    }
    return fp;
}

void makeDir (char *dirName, mode_t dirMode)
{
    if (mkdir(dirName, dirMode) == -1 && errno != EEXIST) {
        sprintf(ErrMsg, "cannot create directory \"%s\"", dirName);
        fatalError();
    }
}

void logWarning (void)
{
    fprintf(WrnFP, "Warning: %s\n", WrnMsg);
    WrnCnt++;
}

void logError (void)
{
    fprintf(ErrFP, "ERROR: %s\n", ErrMsg);
    ErrCnt++;
}

void fatalError (void)
{
    fprintf(stderr, "ERROR: %s\n", ErrMsg);
    exit(1);
}
