Sat Sep 18 15:22:54 CDT 2021
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
$DATA ../../../../data/spa/D/dat46.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12556.9101508846        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9191E+02  2.2136E+02 -1.8610E+01  9.4523E+01  9.9580E+01 -1.7236E+03 -7.2635E+02 -4.8703E+01 -1.1411E+03 -4.2016E+02
            -2.4214E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -580.900058886921        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3584E+00  1.2043E+00  9.9482E-01  1.6878E+00  1.1836E+00  1.7915E+00  1.2335E+00  9.8392E-01  1.2428E+00  1.0646E+00
             1.4641E+01
 PARAMETER:  4.0634E-01  2.8589E-01  9.4811E-02  6.2343E-01  2.6852E-01  6.8304E-01  3.0982E-01  8.3788E-02  3.1740E-01  1.6260E-01
             2.7838E+00
 GRADIENT:   1.8116E+01  2.4006E+01 -2.6122E+00  4.2607E+01 -9.9145E+00  3.4579E+01  2.8231E-01  3.9283E+00  9.1536E+00  3.1661E+00
             1.2171E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -595.407729157886        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2776E+00  7.8765E-01  1.3089E+00  1.9712E+00  1.8532E+00  1.6648E+00  3.4249E+00  4.8436E-01  1.2961E+00  3.6005E+00
             1.2900E+01
 PARAMETER:  3.4499E-01 -1.3870E-01  3.6915E-01  7.7862E-01  7.1692E-01  6.0968E-01  1.3311E+00 -6.2493E-01  3.5935E-01  1.3811E+00
             2.6573E+00
 GRADIENT:   1.8599E+01  2.1045E+01 -9.3725E+00  5.5815E+01 -1.4651E+01  5.4706E-02  1.1231E+01  5.4115E-01  1.9972E+01  1.4241E+01
             7.3150E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -628.991877836427        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0653E+00  2.5329E-01  1.7181E+00  1.6190E+00  3.9006E+00  1.4611E+00  3.0772E+00  3.0186E-01  5.0703E-01  3.5902E+00
             1.0942E+01
 PARAMETER:  1.6325E-01 -1.2732E+00  6.4122E-01  5.8178E-01  1.4611E+00  4.7920E-01  1.2240E+00 -1.0978E+00 -5.7918E-01  1.3782E+00
             2.4926E+00
 GRADIENT:   1.9295E+01  1.1785E-01  1.4498E+01 -3.8005E+01 -3.7688E+00  1.6923E+01  1.8949E+00 -1.9537E-02 -7.9460E-01 -1.2335E+00
            -4.2165E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -669.785424632030        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  8.2714E-01  1.2066E-01  1.8331E-01  1.1177E+00  9.4027E+00  1.1742E+00  1.2102E+00  1.0000E-02  1.0018E-01  2.1145E+00
             1.0823E+01
 PARAMETER: -8.9778E-02 -2.0148E+00 -1.5966E+00  2.1131E-01  2.3410E+00  2.6061E-01  2.9080E-01 -8.3623E+00 -2.2008E+00  8.4880E-01
             2.4817E+00
 GRADIENT:  -1.8638E+01  5.2564E+00 -3.7277E+01  1.5271E+02 -1.8792E+00 -9.3722E+01  2.7362E+00  0.0000E+00  3.2087E-01  2.4671E+00
            -2.1958E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -702.049846710660        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  6.0265E-01  4.5850E-02  6.2093E-02  5.7939E-01  2.8651E+01  1.4301E+00  1.7882E-01  1.0000E-02  2.6632E-02  1.7404E+00
             9.8982E+00
 PARAMETER: -4.0642E-01 -2.9824E+00 -2.6791E+00 -4.4578E-01  3.4552E+00  4.5772E-01 -1.6214E+00 -1.4862E+01 -3.5257E+00  6.5414E-01
             2.3924E+00
 GRADIENT:   1.7507E+01  1.9907E+00 -9.3706E+00  1.5588E+01 -4.3245E-02  2.4343E+00  1.0573E-01  0.0000E+00  1.4675E-02  1.2942E-02
            -3.1408E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -703.903082440274        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      467
 NPARAMETR:  4.7891E-01  2.2134E-02  3.3906E-02  3.7986E-01  4.7555E+01  1.3647E+00  6.3241E-02  1.0000E-02  1.0000E-02  1.6037E+00
             1.0079E+01
 PARAMETER: -6.3625E-01 -3.7106E+00 -3.2842E+00 -8.6795E-01  3.9619E+00  4.1095E-01 -2.6608E+00 -1.8621E+01 -4.6998E+00  5.7232E-01
             2.4105E+00
 GRADIENT:  -2.6319E+00  1.7164E+00 -1.4323E+01  1.8106E+01 -2.7478E-02 -2.5894E+00  1.7114E-03  0.0000E+00  0.0000E+00  4.3753E-04
            -9.7141E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -704.254727190303        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  4.9386E-01  1.9932E-02  3.6759E-02  3.9909E-01  3.4947E+01  1.3765E+00  7.0891E-02  1.0000E-02  1.0000E-02  1.3277E+00
             1.0261E+01
 PARAMETER: -6.0551E-01 -3.8155E+00 -3.2034E+00 -8.1856E-01  3.6538E+00  4.1953E-01 -2.5466E+00 -1.8472E+01 -4.6874E+00  3.8343E-01
             2.4284E+00
 GRADIENT:  -1.5652E+00  6.2681E-01 -2.5301E+00  3.7727E+00 -1.3086E-02 -3.7219E-01  2.9148E-04  0.0000E+00  0.0000E+00  5.2426E-05
             5.6233E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -704.361959079867        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  4.8653E-01  1.3238E-02  3.5622E-02  3.8935E-01  2.4551E+01  1.3791E+00  7.4276E-02  1.0000E-02  1.0000E-02  1.0051E+00
             1.0186E+01
 PARAMETER: -6.2046E-01 -4.2247E+00 -3.2348E+00 -8.4328E-01  3.3007E+00  4.2145E-01 -2.5000E+00 -1.9006E+01 -4.9139E+00  1.0509E-01
             2.4210E+00
 GRADIENT:  -5.5339E-02 -4.9433E-03  5.0864E-02 -5.4573E-02  6.4874E-03  8.8376E-03  4.9032E-06  0.0000E+00  0.0000E+00 -5.9670E-05
             2.5670E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -704.361964187790        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  4.8680E-01  1.3208E-02  3.5673E-02  3.8976E-01  2.4402E+01  1.3792E+00  7.4521E-02  1.0000E-02  1.0000E-02  1.0017E+00
             1.0186E+01
 PARAMETER: -6.1991E-01 -4.2269E+00 -3.2333E+00 -8.4222E-01  3.2947E+00  4.2152E-01 -2.4967E+00 -1.9001E+01 -4.9130E+00  1.0165E-01
             2.4210E+00
 GRADIENT:  -1.9300E-02 -6.3057E-03  1.8217E-02 -5.1758E-03  6.7106E-03  2.7845E-03  4.8427E-06  0.0000E+00  0.0000E+00 -6.0240E-05
             2.1616E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -704.361965199722        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  4.8691E-01  1.3180E-02  3.5697E-02  3.8995E-01  2.4281E+01  1.3793E+00  7.4627E-02  1.0000E-02  1.0000E-02  9.9851E-01
             1.0186E+01
 PARAMETER: -6.1967E-01 -4.2291E+00 -3.2327E+00 -8.4174E-01  3.2897E+00  4.2157E-01 -2.4953E+00 -1.8999E+01 -4.9127E+00  9.8507E-02
             2.4210E+00
 GRADIENT:  -2.5968E-03 -7.2729E-03  4.4797E-03  1.5692E-02  6.8465E-03  3.2111E-03  4.7805E-06  0.0000E+00  0.0000E+00 -6.0600E-05
             1.4341E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -704.361969485244        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  4.8716E-01  1.3115E-02  3.5745E-02  3.9033E-01  2.3943E+01  1.3794E+00  7.4769E-02  1.0000E-02  1.0000E-02  9.8953E-01
             1.0186E+01
 PARAMETER: -6.1917E-01 -4.2340E+00 -3.2313E+00 -8.4075E-01  3.2757E+00  4.2166E-01 -2.4934E+00 -1.8996E+01 -4.9110E+00  8.9478E-02
             2.4210E+00
 GRADIENT:   3.1183E-02 -9.2691E-03 -2.3370E-02  5.8148E-02  7.1518E-03  4.1569E-03  4.6606E-06  0.0000E+00  0.0000E+00 -6.1486E-05
            -4.1233E-04

0ITERATION NO.:   60    OBJECTIVE VALUE:  -704.362141916631        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1529
 NPARAMETR:  4.8842E-01  1.2631E-02  3.5985E-02  3.9230E-01  2.0390E+01  1.3802E+00  7.3769E-02  1.0000E-02  1.0000E-02  8.8923E-01
             1.0187E+01
 PARAMETER: -6.1657E-01 -4.2716E+00 -3.2246E+00 -8.3573E-01  3.1150E+00  4.2223E-01 -2.5068E+00 -1.8974E+01 -4.8771E+00 -1.7397E-02
             2.4211E+00
 GRADIENT:   1.7857E-01 -1.9056E-02 -1.4474E-01  2.5380E-01  9.4553E-03  7.4660E-03  3.7624E-06  0.0000E+00  0.0000E+00 -6.9921E-05
            -6.3972E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -704.493026870228        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1709
 NPARAMETR:  4.8891E-01  1.0000E-02  3.5443E-02  3.9190E-01  6.7399E+00  1.3855E+00  5.8708E-02  1.0000E-02  1.0000E-02  4.0356E-01
             1.0180E+01
 PARAMETER: -6.1557E-01 -4.5609E+00 -3.2398E+00 -8.3674E-01  2.0080E+00  4.2608E-01 -2.7352E+00 -1.9117E+01 -4.6861E+00 -8.0742E-01
             2.4204E+00
 GRADIENT:   5.3312E-02  0.0000E+00 -2.1030E-01  1.2808E-01  2.4650E-02  3.1445E-02  1.9422E-05  0.0000E+00  0.0000E+00  1.0625E-01
             3.8871E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -704.563276392708        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1888
 NPARAMETR:  4.9666E-01  1.0000E-02  3.6710E-02  4.0396E-01  6.2770E+00  1.3892E+00  9.8100E-02  1.0000E-02  1.0000E-02  1.1591E-01
             1.0176E+01
 PARAMETER: -5.9984E-01 -5.6067E+00 -3.2047E+00 -8.0643E-01  1.9369E+00  4.2871E-01 -2.2218E+00 -2.0464E+01 -5.4728E+00 -2.0550E+00
             2.4200E+00
 GRADIENT:   9.7691E-01  0.0000E+00 -1.3459E+00  1.8494E+00 -4.1260E-03 -3.9845E-01  9.8195E-05  0.0000E+00  0.0000E+00  2.3133E-02
            -8.2426E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -704.580170525975        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2063
 NPARAMETR:  4.9185E-01  1.0000E-02  3.5809E-02  3.9659E-01  6.3038E+00  1.3886E+00  2.3187E-01  1.0000E-02  1.0000E-02  1.6781E-02
             1.0178E+01
 PARAMETER: -6.0958E-01 -7.6684E+00 -3.2295E+00 -8.2484E-01  1.9412E+00  4.2831E-01 -1.3616E+00 -2.3505E+01 -7.1860E+00 -3.9875E+00
             2.4202E+00
 GRADIENT:   7.6030E-02  0.0000E+00 -5.6503E-02  7.5898E-02 -1.4949E-02 -1.7427E-02  5.9630E-04  0.0000E+00  0.0000E+00  4.8821E-04
            -4.6522E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -704.580423873227        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2242            RESET HESSIAN, TYPE II
 NPARAMETR:  4.9135E-01  1.0000E-02  3.5715E-02  3.9584E-01  6.3134E+00  1.3884E+00  1.9023E-01  1.0000E-02  1.0000E-02  1.3546E-02
             1.0178E+01
 PARAMETER: -6.1060E-01 -7.8907E+00 -3.2322E+00 -8.2673E-01  1.9427E+00  4.2815E-01 -1.5595E+00 -2.3832E+01 -7.3709E+00 -4.2017E+00
             2.4202E+00
 GRADIENT:   3.9260E+00  0.0000E+00  3.9112E+00  1.8254E+00  8.0401E-02  4.2462E-01  4.0269E-04  0.0000E+00  0.0000E+00  3.1635E-04
             1.9065E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -704.580567433308        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2420
 NPARAMETR:  4.9134E-01  1.0000E-02  3.5718E-02  3.9583E-01  6.3135E+00  1.3885E+00  1.1344E-01  1.0000E-02  1.0000E-02  1.3317E-02
             1.0178E+01
 PARAMETER: -6.1061E-01 -7.8907E+00 -3.2321E+00 -8.2677E-01  1.9427E+00  4.2821E-01 -2.0765E+00 -2.3832E+01 -7.3709E+00 -4.2187E+00
             2.4202E+00
 GRADIENT:   1.4544E-02  0.0000E+00 -2.5107E-02  3.0220E-02  1.0242E-02  5.8837E-03  1.4230E-04  0.0000E+00  0.0000E+00  3.0302E-04
            -5.4953E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -704.580701515972        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2596
 NPARAMETR:  4.9140E-01  1.0000E-02  3.5727E-02  3.9591E-01  6.3083E+00  1.3885E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.0178E+01
 PARAMETER: -6.1049E-01 -7.8907E+00 -3.2319E+00 -8.2656E-01  1.9419E+00  4.2821E-01 -1.1683E+01 -2.3832E+01 -7.3709E+00 -4.5599E+00
             2.4203E+00
 GRADIENT:   5.6360E-03  0.0000E+00  3.0620E-02 -4.4720E-02 -1.3224E-02 -1.4754E-03  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             8.0071E-03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -704.580703913131        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2771
 NPARAMETR:  4.9139E-01  1.0000E-02  3.5726E-02  3.9590E-01  6.3110E+00  1.3885E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.0178E+01
 PARAMETER: -6.1052E-01 -7.8907E+00 -3.2319E+00 -8.2658E-01  1.9423E+00  4.2821E-01 -1.1093E+01 -2.3832E+01 -7.3709E+00 -4.5388E+00
             2.4202E+00
 GRADIENT:  -8.5173E-06  0.0000E+00 -3.9321E-04  4.2019E-03  2.6568E-04  1.9088E-04  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -2.7597E-05

0ITERATION NO.:   98    OBJECTIVE VALUE:  -704.580704552800        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     2879
 NPARAMETR:  4.9139E-01  1.0000E-02  3.5724E-02  3.9589E-01  6.3112E+00  1.3885E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.0178E+01
 PARAMETER: -6.1053E-01 -7.8907E+00 -3.2319E+00 -8.2662E-01  1.9423E+00  4.2820E-01 -1.0250E+01 -2.3832E+01 -7.3709E+00 -4.5087E+00
             2.4202E+00
 GRADIENT:  -3.3652E-03  0.0000E+00  2.2858E-02 -7.0875E-03 -2.1872E-04 -1.1399E-03  0.0000E+00  0.0000E+00  0.0000E+00  1.2337E-05
             1.7010E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2879
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6780E-03  2.2792E-06  1.1444E-04 -2.4811E-04  3.9782E-06
 SE:             2.8591E-02  2.7216E-06  2.7207E-04  3.6441E-04  3.5493E-05
 N:                     100         100         100         100         100

 P VAL.:         9.2538E-01  4.0234E-01  6.7402E-01  4.9595E-01  9.1076E-01

 ETASHRINKSD(%)  4.2163E+00  9.9991E+01  9.9089E+01  9.8779E+01  9.9881E+01
 ETASHRINKVR(%)  8.2549E+00  1.0000E+02  9.9992E+01  9.9985E+01  1.0000E+02
 EBVSHRINKSD(%)  4.3310E+00  9.9988E+01  9.9078E+01  9.8762E+01  9.9851E+01
 EBVSHRINKVR(%)  8.4744E+00  1.0000E+02  9.9991E+01  9.9985E+01  1.0000E+02
 RELATIVEINF(%)  5.2364E+00  9.5548E-08  5.6804E-05  1.0490E-04  1.7199E-05
 EPSSHRINKSD(%)  7.1357E+00
 EPSSHRINKVR(%)  1.3762E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -704.58070455280017     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       30.570122010938007     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.75
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -704.581       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.91E-01  1.00E-02  3.57E-02  3.96E-01  6.31E+00  1.39E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.02E+01
 


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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.16E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -7.66E+03  0.00E+00  4.37E+05
 
 TH 4
+       -3.48E+02  0.00E+00 -5.16E+04  6.93E+03
 
 TH 5
+        5.51E+00  0.00E+00 -2.96E+02  3.58E+01  4.93E-01
 
 TH 6
+        6.24E+00  0.00E+00  3.75E+02 -7.34E+01  2.37E-01  8.45E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -9.86E-01  0.00E+00  7.99E-01  7.09E-01 -3.75E-03 -2.95E-01  0.00E+00  0.00E+00  0.00E+00  6.23E+01
 
 TH11
+       -2.41E+01  0.00E+00  2.85E+02 -2.38E+01 -2.19E-01  1.72E+00  0.00E+00  0.00E+00  0.00E+00  9.83E-03  4.29E+00
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       40.792
Stop Time:
Sat Sep 18 15:23:36 CDT 2021
