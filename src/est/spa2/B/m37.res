Thu Sep 30 04:11:31 CDT 2021
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
$DATA ../../../../data/spa2/B/dat37.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      600
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2142.61691605185        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3299E+02  1.3189E+01  1.3158E+02 -7.4420E+00 -1.0914E+00  5.7564E+01 -2.2859E+01 -5.0123E+02 -1.0822E+02 -5.7899E+00
            -2.3009E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2362.41675842483        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.6963E-01  1.0302E+00  1.0974E+00  9.6532E-01  1.0458E+00  1.0936E+00  1.0757E+00  2.1674E+00  8.1774E-01  9.2854E-01
             9.5089E-01
 PARAMETER:  6.9158E-02  1.2974E-01  1.9294E-01  6.4704E-02  1.4475E-01  1.8949E-01  1.7299E-01  8.7351E-01 -1.0121E-01  2.5858E-02
             4.9646E-02
 GRADIENT:   4.3281E+02  3.3518E+01  2.8106E+01 -2.7538E+01 -3.9784E+01  1.6064E+02 -2.5136E+01 -5.3400E+01 -2.0950E+01 -2.7077E-01
            -4.8021E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2365.74298715908        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:      282
 NPARAMETR:  9.6887E-01  1.0391E+00  1.1082E+00  9.6592E-01  1.0581E+00  1.0976E+00  1.1327E+00  2.1754E+00  8.1553E-01  9.3360E-01
             9.4926E-01
 PARAMETER:  6.8379E-02  1.3835E-01  2.0276E-01  6.5328E-02  1.5652E-01  1.9317E-01  2.2463E-01  8.7723E-01 -1.0391E-01  3.1293E-02
             4.7930E-02
 GRADIENT:  -6.5306E+01 -5.1134E+01  1.8117E+01 -9.8499E+01 -6.9286E+01  3.6759E+01 -3.3591E+01 -1.2261E+02 -2.0226E+01 -1.7956E+00
            -5.1078E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2367.81992556661        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.6914E-01  1.0387E+00  1.1018E+00  9.6881E-01  1.0577E+00  9.9100E-01  1.1320E+00  2.1704E+00  8.1577E-01  9.4571E-01
             9.5218E-01
 PARAMETER:  6.8655E-02  1.3797E-01  1.9691E-01  6.8314E-02  1.5609E-01  9.0956E-02  2.2402E-01  8.7493E-01 -1.0363E-01  4.4186E-02
             5.0997E-02
 GRADIENT:   4.1593E+02  5.1523E+01  2.5651E+01 -4.8296E+00 -2.9075E+01  6.5196E+01 -9.5152E+00 -5.4434E+01 -1.3888E+01  2.1536E+00
            -4.6253E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2369.10086141696        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  9.6914E-01  1.0387E+00  1.0905E+00  9.7763E-01  1.0577E+00  9.8948E-01  1.1320E+00  2.1705E+00  8.1577E-01  9.4320E-01
             9.6104E-01
 PARAMETER:  6.8654E-02  1.3798E-01  1.8659E-01  7.7377E-02  1.5609E-01  8.9421E-02  2.2402E-01  8.7498E-01 -1.0363E-01  4.1527E-02
             6.0261E-02
 GRADIENT:  -8.0587E+01 -4.2765E+01  1.2002E+01 -7.1873E+01 -5.9397E+01 -1.1956E+00 -3.2871E+01 -1.2147E+02 -1.8711E+01  3.1372E-01
            -3.6956E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2370.21210388133        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      769
 NPARAMETR:  9.6914E-01  1.0387E+00  1.0839E+00  9.8899E-01  1.0577E+00  9.9061E-01  1.1320E+00  2.1709E+00  8.1577E-01  9.3747E-01
             9.7143E-01
 PARAMETER:  6.8654E-02  1.3798E-01  1.8052E-01  8.8932E-02  1.5610E-01  9.0562E-02  2.2402E-01  8.7514E-01 -1.0363E-01  3.5433E-02
             7.1018E-02
 GRADIENT:  -8.0713E+01 -3.2882E+01  8.0967E+00 -4.8550E+01 -5.2724E+01 -7.4633E-01 -3.2371E+01 -1.2082E+02 -1.7892E+01  2.1371E-01
            -2.5051E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2374.54911911459        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  9.9458E-01  1.0555E+00  1.0823E+00  1.0090E+00  1.0572E+00  9.6996E-01  1.3293E+00  2.1656E+00  9.2175E-01  9.2841E-01
             9.9514E-01
 PARAMETER:  9.4566E-02  1.5403E-01  1.7912E-01  1.0897E-01  1.5564E-01  6.9499E-02  3.8466E-01  8.7271E-01  1.8515E-02  2.5713E-02
             9.5130E-02
 GRADIENT:   4.2964E+02  1.0843E+02  1.3100E+01  9.0347E+01 -1.3272E+01  4.9054E+01  4.4774E+01 -5.5529E+01  1.6425E+01  3.8136E+00
             4.0874E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2375.29313237595        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  1.0057E+00  1.0329E+00  1.0823E+00  1.0171E+00  1.0572E+00  9.8975E-01  1.3424E+00  2.1663E+00  8.4149E-01  9.1315E-01
             9.9351E-01
 PARAMETER:  1.0572E-01  1.3237E-01  1.7910E-01  1.1694E-01  1.5566E-01  8.9698E-02  3.9443E-01  8.7303E-01 -7.2585E-02  9.1447E-03
             9.3484E-02
 GRADIENT:   2.6752E+00 -1.0945E+00  8.0740E-01 -1.6228E+00 -3.7936E+01  1.7997E+00  6.2133E-01 -1.2080E+02 -2.2471E-01 -8.0112E-02
            -9.5577E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2375.36397634141        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0067E+00  1.0366E+00  1.0824E+00  1.0175E+00  1.0571E+00  9.8364E-01  1.3278E+00  2.1666E+00  8.5055E-01  9.1301E-01
             9.9361E-01
 PARAMETER:  1.0702E-01  1.3591E-01  1.7912E-01  1.1742E-01  1.5567E-01  8.2920E-02  3.8328E-01  8.7366E-01 -6.0874E-02  8.8638E-03
             9.3597E-02
 GRADIENT:   3.0604E+00 -2.1748E+04 -1.6487E+04  1.2584E+04  1.8946E+04 -9.5580E-01 -3.8584E+03  3.2153E+03  5.3293E-01 -6.5587E-02
             7.0953E-02
 NUMSIGDIG:         1.6         2.3         2.3         2.3         2.3         1.3         2.3         2.3         1.1         1.9
                    3.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1254
 NO. OF SIG. DIGITS IN FINAL EST.:  1.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.2757E-04 -9.2575E-03 -3.9048E-02  3.6037E-03 -2.9730E-02
 SE:             2.9960E-02  2.3528E-02  3.2030E-02  2.3659E-02  2.0201E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7530E-01  6.9397E-01  2.2280E-01  8.7894E-01  1.4111E-01

 ETASHRINKSD(%)  1.0000E-10  2.1179E+01  1.0000E-10  2.0739E+01  3.2323E+01
 ETASHRINKVR(%)  1.0000E-10  3.7872E+01  1.0000E-10  3.7176E+01  5.4198E+01
 EBVSHRINKSD(%)  3.3279E-01  2.1126E+01  2.5544E+01  2.2648E+01  3.1566E+01
 EBVSHRINKVR(%)  6.6446E-01  3.7789E+01  4.4564E+01  4.0167E+01  5.3168E+01
 RELATIVEINF(%)  9.9325E+01  2.2167E+01  3.6718E+01  2.2990E+01  1.9766E+01
 EPSSHRINKSD(%)  3.2143E+01
 EPSSHRINKVR(%)  5.3954E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2375.3639763414135     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1272.6377364958064     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2375.364       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.04E+00  1.08E+00  1.02E+00  1.06E+00  9.83E-01  1.33E+00  2.17E+00  8.51E-01  9.13E-01  9.94E-01
 


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
+        1.13E+03
 
 TH 2
+       -4.87E+06  3.72E+06
 
 TH 3
+       -1.07E+02  1.03E+03  1.96E+06
 
 TH 4
+        1.73E+02 -1.28E+03  1.20E+03  5.17E+06
 
 TH 5
+        1.27E+02 -1.34E+03 -2.31E+06 -1.33E+03  2.73E+06
 
 TH 6
+       -3.04E+00 -2.35E+02 -1.70E+02  2.75E+02  1.99E+02  2.05E+02
 
 TH 7
+       -4.07E+01  4.06E+02  7.49E+05  4.81E+02 -8.83E+05 -6.53E+01  2.86E+05
 
 TH 8
+        1.13E+01 -1.01E+02 -1.98E+05 -1.32E+02  2.34E+05  1.74E+01 -7.56E+04  2.01E+04
 
 TH 9
+        7.41E-01 -6.16E+06  5.70E+02 -9.19E+02 -6.49E+02 -2.42E-01  2.42E+02 -4.89E+01  8.92E+01
 
 TH10
+        1.89E-01 -4.54E+01 -1.17E+01  4.64E+01 -1.15E+01 -1.16E-01 -5.11E+00  6.73E+00  8.70E+00  7.62E+01
 
 TH11
+       -7.14E+00 -5.28E+06  4.19E+02 -6.79E+02 -5.00E+02  2.33E+00  1.62E+02 -3.93E+01  1.21E+01  1.19E+01  5.27E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.079
Stop Time:
Thu Sep 30 04:12:09 CDT 2021
