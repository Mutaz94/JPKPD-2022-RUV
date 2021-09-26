Sat Sep 25 12:52:48 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1695.60969266850        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.6842E+01 -9.6299E+01 -2.8602E+01 -8.8847E+01  7.4055E+01  1.2514E+01  5.1566E+00  7.8983E-01  2.6072E+01 -8.0459E+00
            -1.7303E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1704.65243886004        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      133
 NPARAMETR:  1.0382E+00  1.0641E+00  9.7418E-01  1.0304E+00  9.5486E-01  9.2990E-01  9.2156E-01  1.0115E+00  8.0321E-01  1.0117E+00
             1.0103E+00
 PARAMETER:  1.3744E-01  1.6210E-01  7.3844E-02  1.2992E-01  5.3804E-02  2.7324E-02  1.8309E-02  1.1143E-01 -1.1914E-01  1.1166E-01
             1.1020E-01
 GRADIENT:   1.2381E+01 -4.1958E+00  4.1906E-01 -1.1904E+01  4.5935E+00 -1.8707E+01 -4.6716E+00 -3.2546E-01 -1.4162E+01  1.3163E+00
            -5.6162E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1704.95924174262        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      310
 NPARAMETR:  1.0368E+00  1.0689E+00  8.3493E-01  1.0203E+00  8.7484E-01  9.5287E-01  9.5280E-01  9.7264E-01  7.9220E-01  8.7742E-01
             9.9772E-01
 PARAMETER:  1.3613E-01  1.6664E-01 -8.0405E-02  1.2006E-01 -3.3710E-02  5.1727E-02  5.1647E-02  7.2261E-02 -1.3294E-01 -3.0770E-02
             9.7717E-02
 GRADIENT:   5.6885E+00 -8.6025E-02  6.2550E-01 -5.6272E+00 -1.8804E+00 -9.0309E+00 -4.7944E+00  2.5833E+00 -1.2592E+01 -2.8818E+00
            -4.7436E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1705.88567131357        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  1.0346E+00  1.0736E+00  7.9530E-01  1.0154E+00  8.6159E-01  9.7542E-01  9.5524E-01  7.7551E-01  8.5253E-01  9.1010E-01
             1.0035E+00
 PARAMETER:  1.3402E-01  1.7100E-01 -1.2904E-01  1.1532E-01 -4.8973E-02  7.5110E-02  5.4205E-02 -1.5424E-01 -5.9547E-02  5.7964E-03
             1.0352E-01
 GRADIENT:  -2.9651E-01 -9.9531E-01 -3.5550E-01 -8.1595E-01  3.3749E-01  2.0710E-01 -4.5406E-01  2.2734E-01  8.7317E-02  7.1081E-02
            -2.1400E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1705.89952504656        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  1.0341E+00  9.9711E-01  8.3387E-01  1.0660E+00  8.4797E-01  9.7503E-01  1.0188E+00  7.7876E-01  8.2112E-01  9.0950E-01
             1.0032E+00
 PARAMETER:  1.3352E-01  9.7104E-02 -8.1681E-02  1.6395E-01 -6.4915E-02  7.4718E-02  1.1861E-01 -1.5005E-01 -9.7080E-02  5.1368E-03
             1.0320E-01
 GRADIENT:  -2.2050E-01  1.5049E+00  2.6097E-01  1.9063E+00 -4.1427E-01  3.0090E-01 -2.1366E-02 -7.7036E-02 -1.2419E-01 -9.9008E-02
            -1.8983E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1705.91290061597        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  1.0330E+00  8.9284E-01  9.2471E-01  1.1358E+00  8.5033E-01  9.7169E-01  1.0942E+00  8.5617E-01  7.9084E-01  9.2858E-01
             1.0036E+00
 PARAMETER:  1.3248E-01 -1.3345E-02  2.1720E-02  2.2737E-01 -6.2136E-02  7.1284E-02  1.9000E-01 -5.5287E-02 -1.3465E-01  2.5905E-02
             1.0358E-01
 GRADIENT:   5.0479E-01  2.2249E+00  1.4359E+00  1.6584E+00 -3.2792E+00 -4.8740E-01  3.9753E-02  6.6632E-02  5.1032E-02  1.6314E-01
             8.7616E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1705.95383226375        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  1.0304E+00  7.6229E-01  1.0723E+00  1.2262E+00  8.7245E-01  9.7081E-01  1.1965E+00  9.7787E-01  7.5689E-01  9.6657E-01
             1.0035E+00
 PARAMETER:  1.2993E-01 -1.7143E-01  1.6984E-01  3.0395E-01 -3.6449E-02  7.0376E-02  2.7941E-01  7.7620E-02 -1.7854E-01  6.5998E-02
             1.0348E-01
 GRADIENT:  -2.1590E-01  3.2233E+00  9.6632E-01  5.5482E+00 -1.8551E+00  3.6413E-02 -4.2022E-02 -5.6205E-03 -4.8843E-01 -2.1966E-01
            -2.7927E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1706.07394836436        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1201
 NPARAMETR:  1.0268E+00  5.5287E-01  1.2358E+00  1.3616E+00  8.7440E-01  9.6922E-01  1.4011E+00  1.0757E+00  7.1478E-01  1.0071E+00
             1.0036E+00
 PARAMETER:  1.2646E-01 -4.9263E-01  3.1169E-01  4.0865E-01 -3.4219E-02  6.8731E-02  4.3728E-01  1.7299E-01 -2.3579E-01  1.0706E-01
             1.0359E-01
 GRADIENT:  -3.8827E-01  2.0925E+00  1.2656E+00  6.4574E+00 -1.6367E+00  7.2603E-01 -4.0452E-01 -9.8949E-01 -5.7552E-01 -3.4074E-02
            -1.4523E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1706.19621632546        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1378
 NPARAMETR:  1.0245E+00  4.1206E-01  1.3921E+00  1.4533E+00  8.9127E-01  9.6412E-01  1.6513E+00  1.2219E+00  6.8699E-01  1.0302E+00
             1.0050E+00
 PARAMETER:  1.2422E-01 -7.8658E-01  4.3078E-01  4.7383E-01 -1.5113E-02  6.3458E-02  6.0158E-01  3.0037E-01 -2.7544E-01  1.2973E-01
             1.0503E-01
 GRADIENT:   4.9080E-01  1.1513E+00  1.5699E+00  3.8108E+00 -1.4402E+00 -3.1302E-01 -9.2138E-02 -7.7292E-01 -3.1579E-01 -3.6898E-01
             4.8551E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1706.27087585363        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1554
 NPARAMETR:  1.0220E+00  2.8266E-01  1.5192E+00  1.5369E+00  8.9936E-01  9.6169E-01  2.0822E+00  1.3558E+00  6.6198E-01  1.0455E+00
             1.0047E+00
 PARAMETER:  1.2174E-01 -1.1635E+00  5.1816E-01  5.2978E-01 -6.0699E-03  6.0940E-02  8.3343E-01  4.0441E-01 -3.1252E-01  1.4453E-01
             1.0466E-01
 GRADIENT:   2.8406E-02  6.4578E-01 -3.1598E-01  3.2870E+00  4.1008E-02 -3.8253E-01  2.6526E-01  4.1767E-02 -3.7849E-01  5.4331E-02
            -8.9692E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1706.31732955879        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1732
 NPARAMETR:  1.0201E+00  1.7889E-01  1.6971E+00  1.6088E+00  9.2152E-01  9.6129E-01  2.5519E+00  1.5234E+00  6.4571E-01  1.0649E+00
             1.0061E+00
 PARAMETER:  1.1990E-01 -1.6210E+00  6.2892E-01  5.7552E-01  1.8268E-02  6.0517E-02  1.0368E+00  5.2094E-01 -3.3740E-01  1.6287E-01
             1.0603E-01
 GRADIENT:  -9.8995E-02  7.5093E-01  1.4135E+00  6.3311E+00 -2.9665E+00  1.6422E-01  3.5242E-03 -1.9194E-01 -7.4627E-01  1.2888E-01
            -1.6087E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1706.37922003754        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1908
 NPARAMETR:  1.0188E+00  1.0512E-01  1.8102E+00  1.6583E+00  9.3619E-01  9.6029E-01  2.9632E+00  1.6331E+00  6.3862E-01  1.0778E+00
             1.0068E+00
 PARAMETER:  1.1859E-01 -2.1527E+00  6.9345E-01  6.0576E-01  3.4060E-02  5.9482E-02  1.1863E+00  5.9047E-01 -3.4844E-01  1.7492E-01
             1.0677E-01
 GRADIENT:  -3.8037E-01  2.4291E-01 -6.2109E-01  5.3470E+00  5.3677E-01  2.1303E-01 -1.2488E-01  5.0815E-02 -3.0797E-02  1.5919E-02
             1.5816E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1706.43442875999        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2084
 NPARAMETR:  1.0179E+00  3.9520E-02  1.9445E+00  1.7030E+00  9.5096E-01  9.5815E-01  3.7875E+00  1.7505E+00  6.2981E-01  1.0877E+00
             1.0080E+00
 PARAMETER:  1.1772E-01 -3.1310E+00  7.6500E-01  6.3240E-01  4.9719E-02  5.7253E-02  1.4317E+00  6.5991E-01 -3.6234E-01  1.8408E-01
             1.0793E-01
 GRADIENT:  -1.2113E-01  9.8105E-02  2.1368E-01  4.1446E+00 -8.6660E-01 -2.6963E-01 -3.3415E-02  5.6840E-02  6.0776E-01  7.5194E-02
             1.3571E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1706.46056003203        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2259
 NPARAMETR:  1.0176E+00  1.5799E-02  2.0030E+00  1.7192E+00  9.5797E-01  9.5878E-01  4.7124E+00  1.7964E+00  6.2395E-01  1.0909E+00
             1.0079E+00
 PARAMETER:  1.1748E-01 -4.0478E+00  7.9465E-01  6.4184E-01  5.7066E-02  5.7908E-02  1.6502E+00  6.8578E-01 -3.7169E-01  1.8701E-01
             1.0782E-01
 GRADIENT:   5.5940E-02  3.9722E-02  6.6975E-01  2.5264E+00 -9.8860E-01  1.3279E-01 -8.0333E-03 -1.7100E-01 -6.5708E-02 -9.0886E-02
            -1.7459E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1706.47046077948        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2434
 NPARAMETR:  1.0175E+00  1.0000E-02  1.9844E+00  1.7213E+00  9.5376E-01  9.5834E-01  5.2759E+00  1.7846E+00  6.2317E-01  1.0895E+00
             1.0077E+00
 PARAMETER:  1.1732E-01 -4.5395E+00  7.8532E-01  6.4307E-01  5.2660E-02  5.7449E-02  1.7631E+00  6.7917E-01 -3.7294E-01  1.8574E-01
             1.0770E-01
 GRADIENT:   5.1342E-02  0.0000E+00 -6.8522E-02  3.2406E-01  6.1086E-02  6.1753E-03 -3.5304E-03 -8.3782E-03 -5.1223E-02  6.6418E-03
            -2.1475E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1706.47627954899        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2619
 NPARAMETR:  1.0175E+00  1.0000E-02  1.9884E+00  1.7212E+00  9.5455E-01  9.5839E-01  1.0320E+01  1.7878E+00  6.2304E-01  1.0898E+00
             1.0078E+00
 PARAMETER:  1.1735E-01 -4.5388E+00  7.8733E-01  6.4305E-01  5.3489E-02  5.7495E-02  2.4340E+00  6.8099E-01 -3.7314E-01  1.8600E-01
             1.0776E-01
 GRADIENT:   2.6702E-02  0.0000E+00  1.9961E-02 -3.8158E-01  2.7671E-02  9.6100E-03 -9.3067E-03  1.4122E-02  3.7231E-01  4.0772E-02
             7.6959E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1706.59605328497        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2802
 NPARAMETR:  1.0173E+00  1.0000E-02  1.9855E+00  1.7180E+00  9.5661E-01  9.5811E-01  1.1591E+01  1.7864E+00  6.1725E-01  1.0760E+00
             1.0084E+00
 PARAMETER:  1.1715E-01 -4.5388E+00  7.8586E-01  6.4119E-01  5.5642E-02  5.7207E-02  2.5502E+00  6.8019E-01 -3.8249E-01  1.7327E-01
             1.0837E-01
 GRADIENT:  -4.8440E-01  0.0000E+00 -3.7978E-01  7.3303E+00  1.7597E+00 -9.1152E-02 -1.0533E+00 -1.9633E-01  3.0747E+00  3.7758E-01
            -2.7625E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1706.61609936594        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2978
 NPARAMETR:  1.0173E+00  1.0000E-02  1.9678E+00  1.7149E+00  9.4994E-01  9.5819E-01  1.1521E+01  1.7777E+00  6.1315E-01  1.0678E+00
             1.0096E+00
 PARAMETER:  1.1715E-01 -4.5388E+00  7.7691E-01  6.3934E-01  4.8640E-02  5.7293E-02  2.5442E+00  6.7529E-01 -3.8915E-01  1.6561E-01
             1.0960E-01
 GRADIENT:  -4.2578E-01  0.0000E+00  2.1852E-01 -1.9695E+00  1.7477E-01 -5.1432E-02  8.4782E-01  7.5801E-02 -8.7492E-01 -5.5354E-02
             2.0735E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1706.62097733534        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3154
 NPARAMETR:  1.0174E+00  1.0000E-02  1.9208E+00  1.7130E+00  9.3927E-01  9.5829E-01  1.1505E+01  1.7394E+00  6.1428E-01  1.0606E+00
             1.0085E+00
 PARAMETER:  1.1724E-01 -4.5388E+00  7.5272E-01  6.3827E-01  3.7346E-02  5.7397E-02  2.5428E+00  6.5356E-01 -3.8731E-01  1.5886E-01
             1.0846E-01
 GRADIENT:  -4.2013E-02  0.0000E+00  4.1336E-03  3.8683E-02 -9.9562E-03 -4.2000E-03 -2.3115E-02  8.9248E-04  1.5518E-02  2.9070E-03
            -2.1116E-03

0ITERATION NO.:   91    OBJECTIVE VALUE:  -1706.62097733534        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     3176
 NPARAMETR:  1.0174E+00  1.0000E-02  1.9208E+00  1.7130E+00  9.3927E-01  9.5829E-01  1.1505E+01  1.7394E+00  6.1428E-01  1.0606E+00
             1.0085E+00
 PARAMETER:  1.1724E-01 -4.5388E+00  7.5272E-01  6.3827E-01  3.7346E-02  5.7397E-02  2.5428E+00  6.5356E-01 -3.8731E-01  1.5886E-01
             1.0846E-01
 GRADIENT:  -4.2013E-02  0.0000E+00  4.1336E-03  3.8683E-02 -9.9562E-03 -4.2000E-03 -2.3115E-02  8.9248E-04  1.5518E-02  2.9070E-03
            -2.1116E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3176
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1390E-06  5.7834E-03 -3.7740E-02 -1.0996E-02 -4.6828E-02
 SE:             2.9836E-02  4.7764E-03  1.9211E-02  2.8654E-02  1.9842E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9981E-01  2.2596E-01  4.9472E-02  7.0116E-01  1.8275E-02

 ETASHRINKSD(%)  4.4639E-02  8.3999E+01  3.5640E+01  4.0042E+00  3.3526E+01
 ETASHRINKVR(%)  8.9258E-02  9.7440E+01  5.8579E+01  7.8481E+00  5.5812E+01
 EBVSHRINKSD(%)  4.4419E-01  8.7836E+01  3.8873E+01  4.1666E+00  2.9447E+01
 EBVSHRINKVR(%)  8.8641E-01  9.8520E+01  6.2634E+01  8.1595E+00  5.0223E+01
 RELATIVEINF(%)  9.9025E+01  9.9603E-01  1.0303E+01  5.6892E+01  1.3637E+01
 EPSSHRINKSD(%)  4.5215E+01
 EPSSHRINKVR(%)  6.9986E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1706.6209773353362     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -971.47015077159801     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.50
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1706.621       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  1.92E+00  1.71E+00  9.39E-01  9.58E-01  1.15E+01  1.74E+00  6.14E-01  1.06E+00  1.01E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.04E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.25E-02  0.00E+00  1.02E-01
 
 TH 4
+        1.99E+00  0.00E+00 -6.22E+00  3.81E+02
 
 TH 5
+        2.38E-01  0.00E+00 -7.43E-01  4.55E+01  5.44E+00
 
 TH 6
+       -2.29E-03  0.00E+00  7.16E-03 -4.38E-01 -5.24E-02  5.05E-04
 
 TH 7
+        4.92E-03  0.00E+00 -1.54E-02  9.41E-01  1.12E-01 -1.08E-03  2.33E-03
 
 TH 8
+       -1.57E-02  0.00E+00  4.91E-02 -3.01E+00 -3.59E-01  3.46E-03 -7.43E-03  2.38E-02
 
 TH 9
+        2.03E+02  0.00E+00 -6.36E+02  3.89E+04  4.65E+03 -4.48E+01  9.62E+01 -3.08E+02  3.98E+06
 
 TH10
+       -1.38E-01  0.00E+00  4.33E-01 -2.65E+01 -3.17E+00  3.05E-02 -6.55E-02  2.09E-01 -2.71E+03  1.85E+00
 
 TH11
+        1.46E-03  0.00E+00 -4.57E-03  2.80E-01  3.34E-02 -3.22E-04  6.91E-04 -2.21E-03  2.86E+01 -1.95E-02  2.06E-04
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.59E+00  0.00E+00  4.21E+01
 
 TH 4
+        1.46E+01  0.00E+00 -9.99E+01  1.15E+05
 
 TH 5
+        3.03E+01  0.00E+00 -1.27E+02  5.01E+02  6.52E+02
 
 TH 6
+       -3.70E+00  0.00E+00 -1.07E+00 -4.77E+00 -1.64E+01  1.72E+02
 
 TH 7
+       -8.91E-01  0.00E+00  2.77E+00  6.51E+01 -2.07E+01  1.30E-01  2.26E+02
 
 TH 8
+        7.67E-02  0.00E+00 -1.30E+01 -4.13E+01 -2.10E+01 -1.36E+00  1.36E+00  2.08E+01
 
 TH 9
+        1.22E+02  0.00E+00 -3.68E+02 -1.45E+03  2.74E+03 -2.63E+01  4.64E+01 -1.79E+02  2.44E+06
 
 TH10
+       -4.43E+00  0.00E+00  1.20E+01 -3.39E+02 -1.64E+02  2.56E+01  1.20E+01  1.56E+01 -1.59E+03  9.79E+01
 
 TH11
+        1.06E+01  0.00E+00 -1.94E+00 -1.26E+01 -1.96E+01 -5.65E+00 -8.67E-02  4.06E+00  2.03E+01  2.02E+01  1.79E+02
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.16E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.69E+00  0.00E+00  3.90E+01
 
 TH 4
+       -2.73E+02  0.00E+00  1.59E+02  6.67E+04
 
 TH 5
+       -4.95E+00  0.00E+00 -1.09E+02 -1.15E+03  5.92E+02
 
 TH 6
+       -1.42E+02  0.00E+00 -1.11E+00 -2.64E+02  9.55E+00  2.13E+02
 
 TH 7
+        8.58E+00  0.00E+00 -6.90E+00 -2.42E+03  3.95E+01  9.02E+00  8.93E+01
 
 TH 8
+        9.18E+00  0.00E+00 -1.41E+01  3.83E+01 -1.11E+01  2.81E-02 -1.43E+00  1.87E+01
 
 TH 9
+       -1.02E+03  0.00E+00  8.50E+02  2.99E+05 -4.82E+03 -1.19E+03 -1.10E+04  1.72E+02  1.36E+06
 
 TH10
+        2.86E+01  0.00E+00 -9.71E-01  2.16E+02 -5.74E+01 -5.12E+00 -8.20E+00  7.51E+00  9.97E+02  4.72E+01
 
 TH11
+        2.72E+01  0.00E+00 -3.17E+00 -1.63E+02 -1.92E+01 -3.82E+01  8.07E+00 -7.25E-01 -9.57E+02  9.44E+00  2.65E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       47.624
Stop Time:
Sat Sep 25 12:53:37 CDT 2021
