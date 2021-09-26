Sat Sep 25 11:03:42 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/SL2/dat32.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m32.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1721.27145204385        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.6281E+01 -6.4608E+01 -3.0445E+00 -7.6442E+01  4.3471E+01  4.1203E+01  2.5926E+00 -1.2808E+01  1.3989E+01  1.0537E+01
             2.6499E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1731.35096298665        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0515E+00  1.0174E+00  8.9649E-01  1.0234E+00  9.0126E-01  8.1605E-01  9.6622E-01  1.1563E+00  8.7513E-01  8.5918E-01
             8.9743E-01
 PARAMETER:  1.5020E-01  1.1720E-01 -9.2674E-03  1.2312E-01 -3.9630E-03 -1.0327E-01  6.5636E-02  2.4522E-01 -3.3382E-02 -5.1775E-02
            -8.2234E-03
 GRADIENT:   1.0873E+02 -1.1011E+01  7.0050E+00 -2.6869E+01 -7.6552E+00 -3.4428E+01 -4.9950E+00 -2.1913E+00 -8.6691E+00  1.1467E+01
            -1.8638E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1733.76998561109        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0492E+00  9.9498E-01  6.2400E-01  1.0228E+00  7.4612E-01  8.3566E-01  1.1702E+00  8.3379E-01  8.4180E-01  6.2245E-01
             9.1060E-01
 PARAMETER:  1.4805E-01  9.4970E-02 -3.7161E-01  1.2256E-01 -1.9287E-01 -7.9529E-02  2.5716E-01 -8.1773E-02 -7.2217E-02 -3.7409E-01
             6.3501E-03
 GRADIENT:   8.7375E+01 -9.9619E-01 -1.7373E+01  1.8111E+01  2.8141E+01 -2.4926E+01  1.0174E+01 -7.2818E-01  1.3608E-01  4.0671E+00
            -1.2269E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1734.86533711764        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  1.0318E+00  1.0338E+00  5.5928E-01  9.8696E-01  7.1240E-01  8.7220E-01  1.0667E+00  7.6512E-01  8.6512E-01  5.5348E-01
             9.3022E-01
 PARAMETER:  1.3126E-01  1.3319E-01 -4.8110E-01  8.6875E-02 -2.3911E-01 -3.6734E-02  1.6455E-01 -1.6772E-01 -4.4891E-02 -4.9154E-01
             2.7668E-02
 GRADIENT:  -3.7893E+01  4.5562E-01 -8.7358E-02 -2.1664E+00 -3.0601E+00 -1.1906E+01  3.1705E-01 -2.3843E+00 -3.1135E-01  4.4835E-01
            -3.8085E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1735.97832583893        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  1.0431E+00  8.8147E-01  8.0340E-01  1.1115E+00  7.8655E-01  8.9177E-01  1.2233E+00  1.0597E+00  8.1038E-01  6.3132E-01
             9.3164E-01
 PARAMETER:  1.4222E-01 -2.6164E-02 -1.1890E-01  2.0571E-01 -1.4010E-01 -1.4542E-02  3.0158E-01  1.5796E-01 -1.1026E-01 -3.5994E-01
             2.9193E-02
 GRADIENT:   7.0714E-01  4.2741E+00  4.6714E+00 -2.4440E+00 -5.8458E+00 -1.2308E+00  4.4246E-01 -8.9075E-02 -5.6371E-01 -2.1495E-01
            -2.7850E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1736.10956309202        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      615             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0423E+00  8.0887E-01  8.1210E-01  1.1589E+00  7.6786E-01  8.9271E-01  1.3090E+00  1.0437E+00  7.8260E-01  6.2063E-01
             9.3993E-01
 PARAMETER:  1.4143E-01 -1.1212E-01 -1.0813E-01  2.4751E-01 -1.6415E-01 -1.3489E-02  3.6929E-01  1.4282E-01 -1.4513E-01 -3.7701E-01
             3.8047E-02
 GRADIENT:   6.6532E+01  5.8273E+00 -3.8374E+00  3.6277E+01  7.3360E+00  4.0422E+00  8.3126E-01  7.9060E-01  2.7162E-01 -4.1945E-02
             7.9933E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1736.12249026981        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  1.0423E+00  8.0887E-01  8.1666E-01  1.1589E+00  7.6786E-01  8.9387E-01  1.3090E+00  1.0437E+00  7.8260E-01  6.2063E-01
             9.3835E-01
 PARAMETER:  1.4143E-01 -1.1212E-01 -1.0253E-01  2.4751E-01 -1.6415E-01 -1.2191E-02  3.6929E-01  1.4282E-01 -1.4513E-01 -3.7701E-01
             3.6371E-02
 GRADIENT:  -1.2500E-01  4.7206E+00  1.4423E-03  8.1278E+00 -7.1524E-01 -1.3303E-03 -5.3530E-01  1.9415E-01 -7.8672E-01 -2.8662E-01
             3.7112E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1736.15034377177        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.0420E+00  8.0254E-01  8.1704E-01  1.1577E+00  7.6747E-01  8.9366E-01  1.3218E+00  1.0337E+00  7.8472E-01  6.2389E-01
             9.3765E-01
 PARAMETER:  1.4114E-01 -1.1997E-01 -1.0207E-01  2.4640E-01 -1.6466E-01 -1.2429E-02  3.7898E-01  1.3317E-01 -1.4242E-01 -3.7177E-01
             3.5619E-02
 GRADIENT:  -6.7830E-01  4.0679E-01 -1.3747E-01 -1.4397E+00  1.0357E+00 -5.7465E-02  8.5374E-02 -1.0347E-01  2.3916E-01 -4.5495E-02
            -1.3506E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1736.17503215472        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0381E+00  7.7194E-01  8.1976E-01  1.1752E+00  7.5722E-01  8.9198E-01  1.3644E+00  1.0279E+00  7.7179E-01  6.1518E-01
             9.3790E-01
 PARAMETER:  1.3744E-01 -1.5885E-01 -9.8738E-02  2.6140E-01 -1.7810E-01 -1.4316E-02  4.1073E-01  1.2757E-01 -1.5904E-01 -3.8584E-01
             3.5885E-02
 GRADIENT:  -1.0228E+01  5.5660E-01  3.6365E-01 -1.7320E+00 -6.5244E-01 -6.9982E-01 -3.0910E-01 -7.6266E-02 -3.3555E-01 -1.6113E-01
            -1.0395E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1736.17611012763        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.0381E+00  7.7116E-01  8.1959E-01  1.1751E+00  7.5728E-01  8.9340E-01  1.3678E+00  1.0279E+00  7.7280E-01  6.1524E-01
             9.3792E-01
 PARAMETER:  1.3740E-01 -1.5986E-01 -9.8955E-02  2.6134E-01 -1.7802E-01 -1.2718E-02  4.1323E-01  1.2753E-01 -1.5773E-01 -3.8575E-01
             3.5910E-02
 GRADIENT:  -1.0257E+01 -1.8081E-02 -1.8203E-01 -2.5525E+00  2.4399E-01 -6.9066E-02 -4.4673E-02  1.3893E-02  4.8433E-02 -1.3090E-01
            -4.2750E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1736.20713803572        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1431
 NPARAMETR:  1.0409E+00  7.6114E-01  8.2390E-01  1.1818E+00  7.5631E-01  8.9262E-01  1.3791E+00  1.0233E+00  7.7073E-01  6.2027E-01
             9.3856E-01
 PARAMETER:  1.4008E-01 -1.7293E-01 -9.3706E-02  2.6706E-01 -1.7930E-01 -1.3594E-02  4.2141E-01  1.2304E-01 -1.6042E-01 -3.7759E-01
             3.6590E-02
 GRADIENT:  -2.5558E+00 -2.5702E-02 -3.6228E-01 -1.4894E+00  2.5829E-01 -3.1772E-01 -2.6798E-01 -1.7418E-01  2.3312E-01  1.9318E-01
             2.4588E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1736.31746043579        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  1.0378E+00  6.1736E-01  8.9124E-01  1.2663E+00  7.4361E-01  8.9064E-01  1.6460E+00  1.0797E+00  7.2436E-01  6.0825E-01
             9.3899E-01
 PARAMETER:  1.3710E-01 -3.8230E-01 -1.5141E-02  3.3609E-01 -1.9624E-01 -1.5810E-02  5.9838E-01  1.7668E-01 -2.2246E-01 -3.9717E-01
             3.7051E-02
 GRADIENT:  -4.5734E+00 -2.5385E+00 -1.3122E+00 -9.4548E+00  1.2185E+00 -1.9565E-01  1.9496E-01  1.1668E+00 -1.3675E-01 -7.7713E-01
             1.9007E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1736.39286126980        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1782
 NPARAMETR:  1.0391E+00  5.4991E-01  9.2904E-01  1.3128E+00  7.4263E-01  8.9029E-01  1.7957E+00  1.0810E+00  7.0834E-01  6.3116E-01
             9.3818E-01
 PARAMETER:  1.3832E-01 -4.9800E-01  2.6394E-02  3.7219E-01 -1.9756E-01 -1.6206E-02  6.8541E-01  1.7791E-01 -2.4484E-01 -3.6019E-01
             3.6189E-02
 GRADIENT:   2.3957E+00 -1.3948E-02 -1.0463E+00  8.4631E-01  3.5959E-01  2.0177E-01  9.6193E-02 -3.8254E-01 -3.7482E-02 -1.9976E-01
            -2.3202E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1736.44156722630        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  1.0354E+00  4.5880E-01  1.0373E+00  1.3749E+00  7.6750E-01  8.8762E-01  2.0311E+00  1.1813E+00  6.8842E-01  6.6933E-01
             9.3836E-01
 PARAMETER:  1.3481E-01 -6.7914E-01  1.3666E-01  4.1837E-01 -1.6461E-01 -1.9210E-02  8.0857E-01  2.6665E-01 -2.7335E-01 -3.0148E-01
             3.6374E-02
 GRADIENT:  -7.3061E-01 -3.4112E-02 -4.8649E-02 -1.3980E-01 -2.8336E-02 -5.6594E-02  2.3055E-03 -4.0396E-02  1.8706E-02  5.7478E-02
            -3.1835E-02

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1736.44170029532        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     2018
 NPARAMETR:  1.0355E+00  4.5521E-01  1.0424E+00  1.3775E+00  7.6873E-01  8.8763E-01  2.0416E+00  1.1872E+00  6.8757E-01  6.7011E-01
             9.3843E-01
 PARAMETER:  1.3490E-01 -6.8701E-01  1.4153E-01  4.2024E-01 -1.6301E-01 -1.9199E-02  8.1371E-01  2.7159E-01 -2.7459E-01 -3.0031E-01
             3.6453E-02
 GRADIENT:  -1.6268E-01 -5.9475E-03 -3.0951E-02  9.1086E-02  4.7505E-02 -1.8107E-02 -4.5417E-03  1.1073E-02 -1.2734E-02  4.7440E-05
            -1.2612E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2018
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1653E-03  2.3236E-02 -3.3586E-02 -2.0322E-02 -2.5301E-02
 SE:             2.9856E-02  1.8174E-02  2.0163E-02  2.4487E-02  1.7354E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6887E-01  2.0106E-01  9.5774E-02  4.0659E-01  1.4485E-01

 ETASHRINKSD(%)  1.0000E-10  3.9115E+01  3.2450E+01  1.7965E+01  4.1863E+01
 ETASHRINKVR(%)  1.0000E-10  6.2930E+01  5.4370E+01  3.2703E+01  6.6201E+01
 EBVSHRINKSD(%)  4.8798E-01  4.2803E+01  3.2877E+01  1.6392E+01  4.0051E+01
 EBVSHRINKVR(%)  9.7359E-01  6.7285E+01  5.4945E+01  3.0098E+01  6.4061E+01
 RELATIVEINF(%)  9.7744E+01  3.5897E+00  6.5994E+00  8.9501E+00  5.0906E+00
 EPSSHRINKSD(%)  4.5811E+01
 EPSSHRINKVR(%)  7.0636E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1736.4417002953235     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1001.2908737315853     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.39
 Elapsed covariance  time in seconds:     6.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1736.442       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  4.55E-01  1.04E+00  1.38E+00  7.69E-01  8.88E-01  2.04E+00  1.19E+00  6.88E-01  6.70E-01  9.38E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.98E-02  6.40E-01  7.00E-01  4.23E-01  1.61E-01  8.00E-02  1.86E+00  7.73E-01  1.59E-01  2.68E-01  7.06E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.89E-04
 
 TH 2
+        7.65E-03  4.09E-01
 
 TH 3
+       -7.80E-03 -4.20E-01  4.90E-01
 
 TH 4
+       -5.02E-03 -2.69E-01  2.84E-01  1.79E-01
 
 TH 5
+       -1.39E-03 -7.93E-02  1.06E-01  5.52E-02  2.60E-02
 
 TH 6
+        1.72E-04  8.05E-03 -7.59E-03 -5.36E-03 -1.33E-03  6.40E-03
 
 TH 7
+       -2.19E-02 -1.18E+00  1.22E+00  7.77E-01  2.32E-01 -2.27E-02  3.45E+00
 
 TH 8
+       -8.14E-03 -4.40E-01  5.25E-01  2.98E-01  1.15E-01 -7.83E-03  1.28E+00  5.97E-01
 
 TH 9
+        1.84E-03  9.51E-02 -1.02E-01 -6.33E-02 -2.04E-02  1.67E-03 -2.76E-01 -1.06E-01  2.52E-02
 
 TH10
+       -2.26E-03 -1.38E-01  1.59E-01  9.36E-02  3.48E-02 -2.60E-03  4.01E-01  1.56E-01 -3.38E-02  7.19E-02
 
 TH11
+       -1.67E-04 -9.92E-03  1.40E-02  7.29E-03  3.35E-03 -6.81E-04  2.86E-02  1.54E-02 -2.72E-03  3.03E-03  4.99E-03
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.98E-02
 
 TH 2
+        4.01E-01  6.40E-01
 
 TH 3
+       -3.74E-01 -9.39E-01  7.00E-01
 
 TH 4
+       -3.98E-01 -9.94E-01  9.57E-01  4.23E-01
 
 TH 5
+       -2.88E-01 -7.68E-01  9.35E-01  8.08E-01  1.61E-01
 
 TH 6
+        7.20E-02  1.57E-01 -1.35E-01 -1.58E-01 -1.03E-01  8.00E-02
 
 TH 7
+       -3.95E-01 -9.92E-01  9.37E-01  9.88E-01  7.75E-01 -1.53E-01  1.86E+00
 
 TH 8
+       -3.53E-01 -8.91E-01  9.70E-01  9.12E-01  9.20E-01 -1.27E-01  8.91E-01  7.73E-01
 
 TH 9
+        3.89E-01  9.36E-01 -9.16E-01 -9.41E-01 -7.95E-01  1.31E-01 -9.36E-01 -8.63E-01  1.59E-01
 
 TH10
+       -2.82E-01 -8.06E-01  8.48E-01  8.25E-01  8.04E-01 -1.21E-01  8.06E-01  7.54E-01 -7.93E-01  2.68E-01
 
 TH11
+       -7.94E-02 -2.20E-01  2.83E-01  2.44E-01  2.94E-01 -1.20E-01  2.18E-01  2.82E-01 -2.42E-01  1.60E-01  7.06E-02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.37E+03
 
 TH 2
+        5.25E+00  5.53E+02
 
 TH 3
+        7.63E+01  1.76E+02  3.94E+02
 
 TH 4
+        5.68E+00  5.02E+02 -4.70E+01  9.35E+02
 
 TH 5
+       -1.41E+02 -4.80E+02 -6.78E+02 -6.96E+01  1.71E+03
 
 TH 6
+       -8.35E+00 -6.03E+00 -1.33E+01  1.63E+01  1.88E+01  1.63E+02
 
 TH 7
+       -4.66E+00  5.26E+01  4.81E+00 -2.80E+00 -9.07E+00 -9.34E-01  1.95E+01
 
 TH 8
+       -1.12E+01 -1.72E+01 -6.30E+01 -1.41E+01  1.65E-01  6.81E-01 -1.86E+00  4.33E+01
 
 TH 9
+       -4.42E+01 -1.87E+01  2.39E-01  5.75E+01  1.12E+02  1.18E+01  1.24E+01 -2.03E+01  3.85E+02
 
 TH10
+       -2.51E+01 -9.81E+00 -3.18E+01 -3.79E+01 -6.43E+01  1.63E+00 -8.62E-01  3.08E+01 -1.13E+01  7.68E+01
 
 TH11
+       -2.06E+01 -5.77E+01 -4.68E+01 -7.29E+01  4.86E+01  1.96E+01  2.14E+00  1.04E+01  1.12E+01  2.49E+01  2.40E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       30.628
Stop Time:
Sat Sep 25 11:04:14 CDT 2021
