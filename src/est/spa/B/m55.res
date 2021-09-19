Sat Sep 18 08:34:42 CDT 2021
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
$DATA ../../../../data/spa/B/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1647.56881432294        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7466E+01 -9.1366E+01 -3.5039E+01 -6.8277E+01  9.9538E+01  5.4864E+01 -1.4070E+01 -6.0017E+00 -8.1936E+00 -8.2621E+00
            -1.8892E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1659.73755499543        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9803E-01  1.0770E+00  9.6280E-01  1.0155E+00  9.3109E-01  8.4087E-01  1.0961E+00  1.0911E+00  1.0139E+00  9.3096E-01
             1.0191E+00
 PARAMETER:  9.8033E-02  1.7420E-01  6.2094E-02  1.1536E-01  2.8604E-02 -7.3316E-02  1.9175E-01  1.8716E-01  1.1383E-01  2.8459E-02
             1.1887E-01
 GRADIENT:   5.4054E+01  1.0359E+01  4.7934E-01  1.5496E+01  1.4031E+00 -4.8186E+00 -1.2090E+00 -1.2567E+00  1.4606E+00  5.2855E-01
            -1.0697E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1659.93673113227        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9559E-01  1.0325E+00  9.3738E-01  1.0406E+00  9.0077E-01  8.5357E-01  1.2212E+00  1.0818E+00  9.6690E-01  8.7945E-01
             1.0513E+00
 PARAMETER:  9.5581E-02  1.3201E-01  3.5334E-02  1.3979E-01 -4.5053E-03 -5.8330E-02  2.9986E-01  1.7864E-01  6.6339E-02 -2.8461E-02
             1.5006E-01
 GRADIENT:   4.3197E+01  9.8151E+00 -3.4509E+00  1.5875E+01  4.5313E+00  1.2286E+00  4.8062E+00  1.5948E+00  4.0992E-01  2.0198E-02
             2.9383E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1660.03336123675        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      254
 NPARAMETR:  9.8746E-01  1.0247E+00  7.9065E-01  1.0270E+00  8.2623E-01  8.5226E-01  1.2205E+00  8.3706E-01  9.6347E-01  8.0976E-01
             1.0450E+00
 PARAMETER:  8.7385E-02  1.2439E-01 -1.3490E-01  1.2662E-01 -9.0886E-02 -5.9862E-02  2.9926E-01 -7.7855E-02  6.2785E-02 -1.1101E-01
             1.4399E-01
 GRADIENT:  -2.0546E+01  1.3279E+00 -2.3849E+00  1.2208E+00  1.1347E+00 -2.0889E+00  9.8034E-01  1.1760E+00 -1.8665E-01  3.1822E-01
             7.0986E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1660.43039080649        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  9.9858E-01  8.5716E-01  6.3110E-01  1.1051E+00  6.7006E-01  8.5817E-01  1.4408E+00  4.6080E-01  8.7967E-01  6.7226E-01
             1.0433E+00
 PARAMETER:  9.8582E-02 -5.4133E-02 -3.6029E-01  1.9994E-01 -3.0039E-01 -5.2948E-02  4.6522E-01 -6.7480E-01 -2.8213E-02 -2.9710E-01
             1.4234E-01
 GRADIENT:   9.4441E+00  6.4410E+00 -2.6098E+00  9.3596E+00  1.0278E+00  2.3797E-01  1.4166E+00  1.3047E-01 -1.5211E+00  5.0623E-01
            -4.0244E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1660.60973182196        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      607
 NPARAMETR:  9.9502E-01  7.6026E-01  6.5600E-01  1.1544E+00  6.5163E-01  8.5731E-01  1.5629E+00  4.7902E-01  8.5524E-01  6.8121E-01
             1.0403E+00
 PARAMETER:  9.5007E-02 -1.7410E-01 -3.2159E-01  2.4356E-01 -3.2829E-01 -5.3952E-02  5.4656E-01 -6.3602E-01 -5.6369E-02 -2.8389E-01
             1.3953E-01
 GRADIENT:   8.4770E-01 -1.1044E+00 -2.0315E-01 -3.8085E+00 -4.0378E-02  1.9321E-01 -3.3038E-01 -3.2601E-01 -1.0200E-01  1.4440E+00
            -1.0238E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.78593490639        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  9.9280E-01  6.4657E-01  7.6389E-01  1.2373E+00  6.7192E-01  8.5534E-01  1.7830E+00  6.7178E-01  8.2262E-01  6.9812E-01
             1.0422E+00
 PARAMETER:  9.2777E-02 -3.3608E-01 -1.6934E-01  3.1293E-01 -2.9762E-01 -5.6260E-02  6.7831E-01 -2.9783E-01 -9.5266E-02 -2.5936E-01
             1.4131E-01
 GRADIENT:   9.1294E-01  7.7681E-01  1.0696E+00  8.0584E-02 -1.2608E+00  1.7532E-01 -2.0367E-01 -1.4230E-02  1.3020E-02 -6.4505E-01
            -4.9824E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.82567974803        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  9.8998E-01  5.5204E-01  8.7128E-01  1.3050E+00  6.9908E-01  8.5303E-01  1.9976E+00  8.0081E-01  7.9766E-01  7.4218E-01
             1.0423E+00
 PARAMETER:  8.9928E-02 -4.9414E-01 -3.7788E-02  3.6618E-01 -2.5800E-01 -5.8959E-02  7.9193E-01 -1.2213E-01 -1.2607E-01 -1.9817E-01
             1.4147E-01
 GRADIENT:  -5.6039E-01  6.6983E-02  3.3806E-01  7.1701E-02  1.2647E-01 -4.0208E-02  6.6311E-04 -2.2506E-01 -2.7770E-02 -1.9465E-01
            -1.1638E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1660.82630286532        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1137
 NPARAMETR:  9.8996E-01  5.4416E-01  8.8103E-01  1.3106E+00  7.0152E-01  8.5296E-01  2.0166E+00  8.1527E-01  7.9571E-01  7.4616E-01
             1.0425E+00
 PARAMETER:  8.9912E-02 -5.0852E-01 -2.6661E-02  3.7051E-01 -2.5450E-01 -5.9042E-02  8.0141E-01 -1.0423E-01 -1.2852E-01 -1.9282E-01
             1.4159E-01
 GRADIENT:   1.1594E-02  7.7069E-04  1.6355E-02 -3.7421E-02  7.5443E-04 -2.5098E-04  1.3852E-03 -2.9658E-03  7.6188E-04 -5.4969E-03
            -9.7878E-04

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1660.82630286532        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1159
 NPARAMETR:  9.8996E-01  5.4416E-01  8.8103E-01  1.3106E+00  7.0152E-01  8.5296E-01  2.0166E+00  8.1527E-01  7.9571E-01  7.4616E-01
             1.0425E+00
 PARAMETER:  8.9912E-02 -5.0852E-01 -2.6661E-02  3.7051E-01 -2.5450E-01 -5.9042E-02  8.0141E-01 -1.0423E-01 -1.2852E-01 -1.9282E-01
             1.4159E-01
 GRADIENT:   1.1594E-02  7.7069E-04  1.6355E-02 -3.7421E-02  7.5443E-04 -2.5098E-04  1.3852E-03 -2.9658E-03  7.6188E-04 -5.4969E-03
            -9.7878E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1159
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.3608E-04  2.1977E-02 -2.9983E-02 -1.9729E-02 -1.2758E-02
 SE:             2.9786E-02  1.9488E-02  1.5382E-02  2.4667E-02  1.9822E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7761E-01  2.5944E-01  5.1270E-02  4.2382E-01  5.1980E-01

 ETASHRINKSD(%)  2.1223E-01  3.4713E+01  4.8468E+01  1.7362E+01  3.3595E+01
 ETASHRINKVR(%)  4.2401E-01  5.7376E+01  7.3444E+01  3.1709E+01  5.5903E+01
 EBVSHRINKSD(%)  6.4369E-01  3.7174E+01  4.9655E+01  1.6417E+01  3.1379E+01
 EBVSHRINKVR(%)  1.2832E+00  6.0529E+01  7.4654E+01  3.0140E+01  5.2912E+01
 RELATIVEINF(%)  9.7757E+01  4.9409E+00  2.9767E+00  1.0887E+01  4.8706E+00
 EPSSHRINKSD(%)  4.4923E+01
 EPSSHRINKVR(%)  6.9666E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.8263028653193     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -925.67547630158117     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.826       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.90E-01  5.44E-01  8.81E-01  1.31E+00  7.02E-01  8.53E-01  2.02E+00  8.15E-01  7.96E-01  7.46E-01  1.04E+00
 


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
+        1.53E+03
 
 TH 2
+       -2.05E+01  4.00E+02
 
 TH 3
+        2.27E+01  1.70E+02  5.72E+02
 
 TH 4
+       -7.38E+00  3.22E+02 -1.75E+02  7.27E+02
 
 TH 5
+        8.41E-01 -4.17E+02 -9.37E+02  1.92E+02  1.92E+03
 
 TH 6
+        1.74E-01 -2.71E+00 -1.30E+00 -3.32E+00 -3.39E+00  2.61E+02
 
 TH 7
+        1.76E+00  3.50E+01 -2.45E+00 -8.08E+00  2.27E+00  1.84E-01  1.58E+01
 
 TH 8
+        4.37E+00 -7.27E+00 -7.41E+01 -5.26E+00  1.69E+01  2.63E+00  1.69E+00  4.29E+01
 
 TH 9
+       -1.06E+00 -2.53E+01 -1.03E+00  9.10E+00 -5.24E+00 -1.11E+00  4.95E+00  1.89E+00  1.82E+02
 
 TH10
+       -3.97E-01  6.12E+00 -2.54E+01 -3.01E+01 -7.78E+01  5.96E-01  5.47E+00  2.38E+01  4.46E+00  8.73E+01
 
 TH11
+       -8.62E+00 -5.28E+00 -1.06E+01 -8.29E+00 -4.58E+00  1.94E+00  1.05E+00  7.44E+00  8.65E+00  1.63E+01  1.94E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.271
Stop Time:
Sat Sep 18 08:35:03 CDT 2021
