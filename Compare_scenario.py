import numpy as np
import math
import copy

class Compare_scenario:
    """
    a simulator class to run a region
    """

    def __init__(self, init_population, beta_0, beta_1, alpha_0, alpha_1,
            b_i, b_min, t_0, resistant_threshold, Ceftriaxone, DrugB, Clinic, Urethritis, CD, Ertapenem,
            Epididymitis, DGI, EpiP, DGIP, EDP, U_E, U_D, U_ED, J_E, J_D, UR, JA, JB, JM, side, mental, risk, hospital,
                 generated_transfer, resistant_probability_A, initial_infection, auto_recover, asy_treat, sym_treat, treat_duration, 
                 fail_duration, rngvalue, newdrug_time, newdrug_time_2, newdrug_time_3, gamma1, gamma2, gamma3, Diagnose, JR, 
                 strategy_flag, resistant_probability_B, Cmin, rp1, rp2, rp3, epsilon):
        """
        state vector
        [S, AI0, AIA, I0, IA, T0, TA, TM, FA]
        Args:
            beta_0 (float):
            beta_1 (float):
        """
        #print('generated_transfer: ', generated_transfer)
        self.strategy_flag = strategy_flag
        self.newdrug_time = newdrug_time
        self.newdrug_time_2 = newdrug_time_2
        self.newdrug_time_3 = newdrug_time_3
        self.rngvalue = rngvalue
        self.auto_recover = auto_recover
        self.asy_treat = asy_treat
        self.sym_treat = sym_treat
        self.treat_duration = treat_duration
        self.fail_duration = fail_duration

        self.M_flag = False
        self.init_transfer = None
        self.rprobA = resistant_probability_A
        self.rprobB = resistant_probability_B
        self.pos_s_0 = initial_infection
        self.stage = 0
        self.judge = 0
        self.state = [[init_population * 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [init_population * 0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [init_population * 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        self.totalstate = [init_population, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.beta_0 = beta_0
        self.beta_1 = beta_1
        self.alpha_0 = alpha_0
        self.alpha_1 = alpha_1
        self.b_i = b_i
        self.b_min = b_min
        self.t_0 = t_0
        self.resistant_threshold = resistant_threshold
        self.CA1 = Ceftriaxone + Urethritis #symptomatic treatment
        self.CA2 = Ceftriaxone + Clinic #asymptomatic treatment
        self.Diagnose = Diagnose
        #print('CA: ', self.CA)
        self.CB1 = DrugB + Urethritis
        self.CB2 = DrugB + Clinic
        self.CS = Epididymitis * EpiP + DGI * DGIP+(Epididymitis + DGI) * EDP
        self.CD = CD
        self.CM = Ertapenem + Clinic
        self.MoriCost = Ertapenem + Clinic
        #self.MoriCost = self.CA
        #self.CM = self.CA
        self.QS = EpiP * U_E * J_E + DGIP * U_D * J_D + EDP * U_ED * J_D
        self.UR = UR
        self.JR = JR
        self.UM = side + mental + risk + hospital
        self.MoriQALY = side + mental + risk + hospital
        #self.MoriQALY = UR
        #self.UM = UR
        self.JA = treat_duration/365
        #self.JA = JA
        #self.JB = JB
        self.JB = treat_duration/365
        self.JM = JM

        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3

        # define transmission parameter for different layers
        
        self.Cmin = Cmin
        self.rp1 = rp1
        self.rp2 = rp2
        self.rp3 = rp3
        self.Clow = self.Cmin * rp1
        self.Cmedium = self.Cmin * rp2
        self.Chigh = self.Cmin * rp3
        self.epsilon = epsilon
        
    






        self.resistant_rate = []
        self.prevalence = []
        self.symptomatic = []
        self.annual_case_rate = []
        self.annual_failtreatment_rate = []
        self.annual_DST_rate = []
        self.annual_case = []
        self.annual_incidence_rate = []
        self.annual_resistance_rateA = []
        self.annual_resistance_rateB = []
        self.annual_resistance_rateAB = []
        self.annual_resistance_rateC = []
        self.annual_resistance_rateBC = []
        self.annual_symptomatic_rate = []
        self.S = []
        self.AI = []
        self.AR = []
        self.I_0 = []
        self.R_0 = []
        self.T_I = []
        self.T_A = []
        self.T_M = []
        self.F_A = []

        self.AI_S = []
        self.AI_T0 = []
        self.AR_S = []
        self.AR_TA = []
        self.I_S = []
        self.I_T0 = []
        self.R_S = []
        self.R_TA = []
        self.T0_S = []
        self.T0_TA = []
        self.TA_F = []
        self.TM_S = []
        self.F_TM = []

        self.Average_CaseRate = []
        self.Average_FailtreatmentRate = []
        self.Average_DSTRate = []
        self.Sum_Cost = []
        self.Sum_QALY = []

        self.Sum_A_cost = []
        self.Sum_B_cost = []
        self.Sum_M_cost = []
        self.Sum_Diagnose_cost = []
        self.Sum_Drugtest_cost = []
        self.Sum_Sequalae_cost = []

        self.Sum_A_qalysloss = []
        self.Sum_B_qalysloss = []
        self.Sum_M_qalysloss = []
        self.Sum_Sequalae_qalysloss = []
        self.Sum_Infection_qalyloss = []




        self.FA_Cost = []
        self.FA_QALY = []
        self.FB_Cost = []
        self.FB_QALY = []
        self.AB_Cost = []
        self.AB_QALY = []
        self.Back_number =[]

        self.cumu_N_0 = []
        self.cumu_I_0 = []
        self.cumu_N_A = []
        self.cumu_I_A = []
        self.cumu_N_B = []
        self.cumu_I_B = []
        self.cumu_N_AB = []
        self.cumu_I_AB = []
        self.cumu_PT_0A = []
        self.cumu_T_0A = []
        self.cumu_PT_AA = []
        self.cumu_T_AA = []
        self.cumu_PT_BA = []
        self.cumu_T_BA = []
        self.cumu_PT_ABA = []
        self.cumu_T_ABA = []
        self.cumu_F_A = []
        self.cumu_F_ABA = []
        self.cumu_PT_0B = []
        self.cumu_T_0B = []
        self.cumu_PT_AB = []
        self.cumu_T_AB = []
        self.cumu_PT_BB = []
        self.cumu_T_BB = []
        self.cumu_PT_ABB = []
        self.cumu_T_ABB = []
        self.cumu_F_BB = []
        self.cumu_F_ABB = []
        self.cumu_finalfail = []
        self.rng = np.random.RandomState(seed=self.rngvalue)
        self.comeyear = self.rng.randint(30, 50)

        self.annul_Cost = []
        self.annul_QALY = []





        '''
        self.transfer_template = [
                [0, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, (618.8, 1710.8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (618.8, 1710.8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (618.8, 1710.8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (618.8, 1710.8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(238, 1643.6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (82242.82885, 5957966.267), 0, 0],
                [0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14), 0],
                [0, 0, 0, 0, 0, 0, 0, (1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1, 14)],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [(1, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        '''

        self.transfer_template = [
                [0, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, self.asy_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.sym_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.asy_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.sym_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.asy_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.sym_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.asy_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.auto_recover, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.sym_treat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.fail_duration, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.fail_duration, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3420104.545, 3420104.545), 0, 0, 0],
                [0, 0, 0, 0, 0, self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.treat_duration, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, self.treat_duration, 0],
                [self.fail_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.fail_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [self.treat_duration, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]



        self.transfer = copy.deepcopy(self.transfer_template)
        self.setup_transfer()
        #previous matrix store
        '''
        if generated_transfer:
            self.transfer = copy.deepcopy(generated_transfer)
        else:
            self.setup_transfer()
            self.init_transfer = copy.deepcopy(self.transfer)
        '''
        #print('self.transfer: ', self.transfer)

    def update_transfer(self):
        """
        update transfer matrix based on transfer template
        for now, we only need to update the first five numbers
        """
        #low activity group
        P_ll = self.epsilon + (1 - self.epsilon) * self.Clow * sum(self.state[0]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        P_lm = (1 - self.epsilon) * self.Cmedium * sum(self.state[1]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        P_lh = (1 - self.epsilon) * self.Chigh * sum(self.state[2]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        # medium activity group
        P_ml = (1 - self.epsilon) * self.Clow * sum(self.state[0]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        P_mm = self.epsilon + (1 - self.epsilon) * self.Cmedium * sum(self.state[1]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        P_mh = (1 - self.epsilon) * self.Chigh * sum(self.state[2]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        # high activity group
        P_hl = (1 - self.epsilon) * self.Clow * sum(self.state[0]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        P_hm = (1 - self.epsilon) * self.Cmedium * sum(self.state[1]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        P_hh = self.epsilon + (1 - self.epsilon) * self.Chigh * sum(self.state[2]) / (self.Clow * sum(self.state[0]) + self.Cmedium * sum(self.state[1]) + self.Chigh * sum(self.state[2]))
        # low activity overall rate
        B_ll = self.beta_0 * self.Clow * P_ll
        B_lm = self.beta_0 * self.Clow * P_lm
        B_lh = self.beta_0 * self.Clow * P_lh
        # medium activity overall rate
        B_ml = self.beta_0 * self.Cmedium * P_ml
        B_mm = self.beta_0 * self.Cmedium * P_mm
        B_mh = self.beta_0 * self.Cmedium * P_mh
        # high activity overall rate
        B_hl = self.beta_0 * self.Chigh * P_hl
        B_hm = self.beta_0 * self.Chigh * P_hm
        B_hh = self.beta_0 * self.Chigh * P_hh
        #define the initial rate spaces
        pos_s_0 = [0, 0, 0]
        pos_s_1 = [0, 0, 0]
        pos_s_2 = [0, 0, 0]
        pos_s_3 = [0, 0, 0]

        pos_s_0[0] = B_ll * (self.state[0][1] + self.state[0][2]) / sum(self.state[0]) + B_lm * (self.state[1][1] + self.state[1][2]) / sum(self.state[1]) + B_lh * (self.state[2][1] + self.state[2][2]) / sum(self.state[2])
        pos_s_0[1] = B_ml * (self.state[0][1] + self.state[0][2]) / sum(self.state[0]) + B_mm * (self.state[1][1] + self.state[1][2]) / sum(self.state[1]) + B_mh * (self.state[2][1] + self.state[2][2]) / sum(self.state[2])
        pos_s_0[2] = B_hl * (self.state[0][1] + self.state[0][2]) / sum(self.state[0]) + B_hm * (self.state[1][1] + self.state[1][2]) / sum(self.state[1]) + B_hh * (self.state[2][1] + self.state[2][2]) / sum(self.state[2])

        pos_s_1[0] = B_ll * self.gamma1 * (self.state[0][3] + self.state[0][4] + self.state[0][17]) / sum(self.state[0]) + B_lm * self.gamma1 * (self.state[1][3] + self.state[1][4] + self.state[1][17]) / sum(self.state[1]) + B_lh * self.gamma1 * (self.state[2][3] + self.state[2][4] + self.state[2][17]) / sum(self.state[2])
        pos_s_1[1] = B_ml * self.gamma1 * (self.state[0][3] + self.state[0][4] + self.state[0][17]) / sum(self.state[0]) + B_mm * self.gamma1 * (self.state[1][3] + self.state[1][4] + self.state[1][17]) / sum(self.state[1]) + B_mh * self.gamma1 * (self.state[2][3] + self.state[2][4] + self.state[2][17]) / sum(self.state[2])
        pos_s_1[2] = B_hl * self.gamma1 * (self.state[0][3] + self.state[0][4] + self.state[0][17]) / sum(self.state[0]) + B_hm * self.gamma1 * (self.state[1][3] + self.state[1][4] + self.state[1][17]) / sum(self.state[1]) + B_hh * self.gamma1 * (self.state[2][3] + self.state[2][4] + self.state[2][17]) / sum(self.state[2])

        pos_s_2[0] = B_ll * self.gamma2 * (self.state[0][5] + self.state[0][6] + self.state[0][27]) / sum(self.state[0]) + B_lm * self.gamma2 * (self.state[1][5] + self.state[1][6] + self.state[1][27]) / sum(self.state[1]) + B_lh * self.gamma2 * (self.state[2][5] + self.state[2][6] + self.state[2][27]) / sum(self.state[2])
        pos_s_2[1] = B_ml * self.gamma2 * (self.state[0][5] + self.state[0][6] + self.state[0][27]) / sum(self.state[0]) + B_mm * self.gamma2 * (self.state[1][5] + self.state[1][6] + self.state[1][27]) / sum(self.state[1]) + B_mh * self.gamma2 * (self.state[2][5] + self.state[2][6] + self.state[2][27]) / sum(self.state[2])
        pos_s_2[2] = B_hl * self.gamma2 * (self.state[0][5] + self.state[0][6] + self.state[0][27]) / sum(self.state[0]) + B_hm * self.gamma2 * (self.state[1][5] + self.state[1][6] + self.state[1][27]) / sum(self.state[1]) + B_hh * self.gamma2 * (self.state[2][5] + self.state[2][6] + self.state[2][27]) / sum(self.state[2])

        pos_s_3[0] = B_ll * self.gamma3 * (self.state[0][7] + self.state[0][8] + self.state[0][18] + self.state[0][28]) / sum(self.state[0]) + B_lm * self.gamma3 * (self.state[1][7] + self.state[1][8] + self.state[1][18] + self.state[1][28]) / sum(self.state[1]) + B_lh * self.gamma3 * (self.state[2][7] + self.state[2][8] + self.state[2][18] + self.state[2][28]) / sum(self.state[2])
        pos_s_3[1] = B_ml * self.gamma3 * (self.state[0][7] + self.state[0][8] + self.state[0][18] + self.state[0][28]) / sum(self.state[0]) + B_mm * self.gamma3 * (self.state[1][7] + self.state[1][8] + self.state[1][18] + self.state[1][28]) / sum(self.state[1]) + B_mh * self.gamma3 * (self.state[2][7] + self.state[2][8] + self.state[2][18] + self.state[2][28]) / sum(self.state[2])
        pos_s_3[2] = B_hl * self.gamma3 * (self.state[0][7] + self.state[0][8] + self.state[0][18] + self.state[0][28]) / sum(self.state[0]) + B_hm * self.gamma3 * (self.state[1][7] + self.state[1][8] + self.state[1][18] + self.state[1][28]) / sum(self.state[1]) + B_hh * self.gamma3 * (self.state[2][7] + self.state[2][8] + self.state[2][18] + self.state[2][28]) / sum(self.state[2])

        #original one: pos_s_0 = self.beta_0 * (self.totalstate[1] + self.totalstate[2]) / sum(self.totalstate) # non-resistance infections
        # self.gamma will be updated during each week old
        '''
        pos_s_1 = self.gamma * self.beta_0 * (self.state[3] + self.state[4] + self.state[17]) / sum(self.state) # A-resistance infections
        pos_s_2 = self.gamma * self.beta_0 * (self.state[5] + self.state[6]) / sum(self.state) # B-resistance infections
        pos_s_3 = self.gamma * self.beta_0 * (self.state[7] + self.state[8] + self.state[18] + self.state[28]) / sum(self.state) # AB-resistance infections
        '''
        '''
        print('pos_s_0: ', pos_s_0)
        print('pos_s_1: ', pos_s_1)
        print('pos_s_2: ', pos_s_2)
        print('pos_s_3: ', pos_s_3)
        '''


         # self.gamma will be updated during each week new
        '''
        pos_s_1 = self.gamma1 * self.beta_0 * (self.totalstate[3] + self.totalstate[4] + self.totalstate[17]) / sum(self.totalstate) # A-resistance infections
        pos_s_2 = self.gamma2 * self.beta_0 * (self.totalstate[5] + self.totalstate[6] + self.totalstate[27]) / sum(self.totalstate) # B-resistance infections
        pos_s_3 = self.gamma3 * self.beta_0 * (self.totalstate[7] + self.totalstate[8] + self.totalstate[18] + self.totalstate[28]) / sum(self.totalstate) # AB-resistance infections
        '''
        # no gamma
        '''
        pos_s_1 = self.beta_0 * (self.state[3] + self.state[4] + self.state[17]) / sum(self.state) # A-resistance infections
        pos_s_2 = self.beta_0 * (self.state[5] + self.state[6] + self.state[27]) / sum(self.state) # B-resistance infections
        pos_s_3 = self.beta_0 * (self.state[7] + self.state[8] + self.state[18] + self.state[28]) / sum(self.state) # AB-resistance infections
        '''
        tmp = [0, 0, 0]
        self.s_0 = [0, 0, 0]
        self.s_a = [0, 0, 0]
        self.s_b = [0, 0, 0]
        self.s_ab = [0, 0, 0]
        for j in range(3):
            tmp[j] = math.e ** -(pos_s_0[j] + pos_s_1[j] + pos_s_2[j] + pos_s_3[j])
            self.s_0[j] = (1 - tmp[j]) * (pos_s_0[j] / (pos_s_0[j] + pos_s_1[j] + pos_s_2[j] + pos_s_3[j]))
            self.s_a[j] = (1 - tmp[j]) * (pos_s_1[j] / (pos_s_0[j] + pos_s_1[j] + pos_s_2[j] + pos_s_3[j]))
            self.s_b[j] = (1 - tmp[j]) * (pos_s_2[j] / (pos_s_0[j] + pos_s_1[j] + pos_s_2[j] + pos_s_3[j]))
            self.s_ab[j] = (1 - tmp[j]) * (pos_s_3[j] / (pos_s_0[j] + pos_s_1[j] + pos_s_2[j] + pos_s_3[j]))

        #print('tmp: ', tmp)
        #self.transfer[0][0] = tmp

        #print('s_0: ', s_0)
        #print('s_a: ', s_a)
        #print('proportion: ', s_a/s_0)




        #need to be removed for different layers
        '''
        self.transfer[0][1] = s_0 * (1 - self.alpha_0)
        self.transfer[0][2] = s_0 * self.alpha_0
        self.transfer[0][3] = s_a * (1 - self.alpha_0)
        self.transfer[0][4] = s_a * self.alpha_0
        self.transfer[0][5] = s_b * (1 - self.alpha_0)
        self.transfer[0][6] = s_b * self.alpha_0
        self.transfer[0][7] = s_ab * (1 - self.alpha_0)
        self.transfer[0][8] = s_ab * self.alpha_0
        self.transfer[0][0] = 1 - sum(self.transfer[0][1:])
        '''


    def setup_transfer(self):
        #rprob = np.random.uniform(low=10**(-6), high=10**(-4))
        rprobA = self.rprobA
        rprobB = self.rprobB
        for i in range(len(self.transfer_template)):
            for j in range(len(self.transfer_template[i])):
                if type(self.transfer_template[i][j]) == type((0, 0)):
                    low = self.transfer_template[i][j][0]
                    high = self.transfer_template[i][j][1]
                    # generate random number for this cell
                    self.transfer[i][j] = np.random.uniform(low=low, high=high)

        #redifine risistance rate according to probability old
        '''
        self.transfer[9][11] = self.transfer[9][0] * (1-rprob) / rprob
        self.transfer[10][12] = self.transfer[10][0] * (1-rprob) / rprob
        self.transfer[13][15] = self.transfer[13][0] * (1-rprob) / rprob
        self.transfer[14][16] = self.transfer[14][0] * (1-rprob) / rprob
        self.transfer[19][23] = self.transfer[19][0] * (1-rprob) / rprob
        self.transfer[20][24] = self.transfer[20][0] * (1-rprob) / rprob
        self.transfer[21][25] = self.transfer[21][0] * (1-rprob) / rprob
        self.transfer[22][26] = self.transfer[22][0] * (1-rprob) / rprob
        '''

        #redifine risistance rate according to probability new
        self.transfer[9][11] = (self.transfer[9][0] * (1-rprobA) / rprobA) * (1 - self.alpha_0)
        self.transfer[9][12] = (self.transfer[9][0] * (1-rprobA) / rprobA) * self.alpha_0
        self.transfer[10][12] = self.transfer[10][0] * (1-rprobA) / rprobA
        self.transfer[13][15] = (self.transfer[13][0] * (1-rprobA) / rprobA) * (1 - self.alpha_0)
        self.transfer[13][16] = (self.transfer[13][0] * (1-rprobA) / rprobA) * self.alpha_0
        self.transfer[14][16] = self.transfer[14][0] * (1-rprobA) / rprobA
        self.transfer[19][23] = (self.transfer[19][0] * (1-rprobB) / rprobB) * (1 - self.alpha_0)
        self.transfer[19][24] = (self.transfer[19][0] * (1-rprobB) / rprobB) * self.alpha_0
        self.transfer[20][24] = self.transfer[20][0] * (1-rprobB) / rprobB
        self.transfer[21][25] = (self.transfer[21][0] * (1-rprobB) / rprobB) * (1 - self.alpha_0)
        self.transfer[21][26] = (self.transfer[21][0] * (1-rprobB) / rprobB) * self.alpha_0
        self.transfer[22][26] = self.transfer[22][0] * (1-rprobB) / rprobB
        if self.strategy_flag == 1:
            self.transfer[27][28] = self.transfer[27][0] * (1-rprobA) / rprobA #add resistance for the final treatment with drug A


        '''
        #keep the recover time same without treatment
        self.transfer[2][0] = self.transfer[1][0]
        self.transfer[3][0] = self.transfer[1][0]
        self.transfer[4][0] = self.transfer[1][0]
        #keep the durtion before the treatment for asymptomatic same
        self.transfer[2][6] = self.transfer[1][5]
        #keep the durtion before the treatment for symptomatic same
        self.transfer[4][6] = self.transfer[3][5]
        '''




        """
        pos_s_0 = np.random.uniform(low=0.0018, high=0.0097)
        pos_s_1 = np.random.uniform(low=0.0001, high=0.000018)
        """
        #pos_s_0 = np.random.uniform(low=0.05, high=0.05)
        #pos_s_0 = np.random.uniform(low=0.01, high=0.1)
        '''
        pos_s_0 = self.pos_s_0
        #pos_s_1 = np.random.uniform(low=0.0001, high=0.001)
        pos_s_1 = np.random.uniform(low=0, high=0)
        pos_s_2 = np.random.uniform(low=0, high=0)
        pos_s_3= np.random.uniform(low=0, high=0)
        '''
        self.s_0 = [0, 0, 0]
        self.s_a = [0, 0, 0]
        self.s_b = [0, 0, 0]
        self.s_ab = [0, 0, 0]
        trans = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]


        #proportion = [0.3, 0.6, 0.1]
        for j in range(3):
            self.s_0[j] = self.pos_s_0 #* proportion[j]
            # pos_s_1 = np.random.uniform(low=0.0001, high=0.001)
            self.s_a[j] = np.random.uniform(low=0, high=0)
            self.s_b[j] = np.random.uniform(low=0, high=0)
            self.s_ab[j] = np.random.uniform(low=0, high=0)
            trans[0][j] = self.s_0[j] * (1 - self.alpha_1)
            trans[1][j] = self.s_0[j] * self.alpha_1
            trans[2][j] = self.s_a[j] * (1 - self.alpha_1)
            trans[3][j] = self.s_a[j] * self.alpha_1
            trans[4][j] = self.s_b[j] * (1 - self.alpha_1)
            trans[5][j] = self.s_b[j] * self.alpha_1
            trans[6][j] = self.s_ab[j] * (1 - self.alpha_1)
            trans[7][j] = self.s_ab[j] * self.alpha_1
        #print('asdasdadasad: ', self.s_0)

        '''
        self.transfer[0][1] = pos_s_0 * (1 - self.alpha_1)
        self.transfer[0][2] = pos_s_0 * self.alpha_1
        self.transfer[0][3] = pos_s_1 * (1 - self.alpha_1)
        self.transfer[0][4] = pos_s_1 * self.alpha_1
        self.transfer[0][5] = pos_s_2 * (1 - self.alpha_1)
        self.transfer[0][6] = pos_s_2 * self.alpha_1
        self.transfer[0][7] = pos_s_3 * (1 - self.alpha_1)
        self.transfer[0][8] = pos_s_3 * self.alpha_1
        '''

        initialproportion = [0.2, 1, 4]  # give different initial prevalence for each group

        for i in range (1,9):
            for j in range(3):
                self.state[j][i] = self.state[j][0] * trans[i-1][j] * initialproportion[j]
                
        
        for j in range(3):
            self.state[j][9] = self.state[j][1]*0.04
            self.state[j][10] = self.state[j][2] * 0.5
            self.state[j][0] = self.state[j][0] - sum(self.state[j][1:])

        #self.state[5] = self.state[0]* pos_s_5
        #self.transfer[0][5] = pos_s_5
        #self.transfer[0][6] = pos_s_6
        #self.transfer[0][7] = pos_s_7
        #self.transfer[0][8] = pos_s_8
        self.transfer[0][0] = 1 - sum(self.transfer[0][1:])

        transfer_copy = self.transfer

        self.AI_S.append(transfer_copy[1][0])
        self.AI_T0.append(transfer_copy[1][5])
        self.AR_S.append(transfer_copy[2][0])
        self.AR_TA.append(transfer_copy[2][6])
        self.I_S.append(transfer_copy[3][0])
        self.I_T0.append(transfer_copy[3][5])
        self.R_S.append(transfer_copy[4][0])
        self.R_TA.append(transfer_copy[4][6])
        self.T0_S.append(transfer_copy[5][0])
        #self.T0_TA.append(transfer_copy[5][6])
        self.T0_TA.append(self.rprobA)
        self.TA_F.append(transfer_copy[6][8])
        self.TM_S.append(transfer_copy[7][0])
        self.F_TM.append(transfer_copy[8][7])


        for i in range(1, len(self.transfer)):
            tmp_pos = [7 / i if i != 0 else 0 for i in self.transfer[i]]
            self.transfer[i] = tmp_pos
            sum_pos = sum(tmp_pos)
            self.transfer[i][i] = math.e ** -sum_pos
            for j in range(len(self.transfer[i])):
                if self.transfer[i][j] != 0 and i != j:
                    self.transfer[i][j] = (1 - self.transfer[i][i]) * (self.transfer[i][j] / sum_pos)

        # self.print_matrix(self.transfer)




    def update_state(self):
        """
        update the state vector based on transfer matrix
        Return:
            T_0_in, T_A_in
        """
        self.init_state = copy.deepcopy(self.state)
        self.state = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        self.totalstate = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        PT0_in = [0, 0, 0]
        T0_in = [0, 0, 0]
        PTA_in = [0, 0, 0]
        TA_in = [0, 0, 0]
        PTB_in = [0, 0, 0]
        TB_in = [0, 0, 0]
        PTAB_in = [0, 0, 0]
        TAB_in = [0, 0, 0]
        TC_in = [0, 0, 0]

        treatment_fail = [0, 0, 0]
        DST_in = [0, 0, 0]


        inci_in = [0, 0, 0]

        Totalcost = [0, 0, 0]
        QALYlost = [0, 0, 0]
        FailAcost = [0, 0, 0]
        FailAQALY = [0, 0, 0]
        FailBcost = [0, 0, 0]
        FailBQALY = [0, 0, 0]
        ABcost = [0, 0, 0]
        ABQALY = [0, 0, 0]
        BeforeswitchAB = [0, 0, 0]
        AfterswitchAB = [0, 0, 0]
        BswitchAB = [0, 0, 0]


        A_cost = [0, 0, 0]
        B_cost = [0, 0, 0]
        M_cost = [0, 0, 0]
        Diagnose_cost = [0, 0, 0]
        Drugtest_cost = [0, 0, 0]
        Sequalae_cost = [0, 0, 0]

        A_qalysloss = [0, 0, 0]
        B_qalysloss = [0, 0, 0]
        M_qalysloss = [0, 0, 0]
        Sequalae_qalysloss = [0, 0, 0]
        Infection_qalyloss = [0, 0, 0]



        N_0 = [0, 0, 0]
        I_0 = [0, 0, 0]
        N_A = [0, 0, 0]
        I_A = [0, 0, 0]
        N_B = [0, 0, 0]
        I_B = [0, 0, 0]
        N_AB = [0, 0, 0]
        I_AB = [0, 0, 0]
        PT_0A = [0, 0, 0]
        T_0A = [0, 0, 0]
        PT_AA = [0, 0, 0]
        T_AA = [0, 0, 0]
        PT_BA = [0, 0, 0]
        T_BA = [0, 0, 0]
        PT_ABA = [0, 0, 0]
        T_ABA = [0, 0, 0]
        F_A = [0, 0, 0]
        F_ABA = [0, 0, 0]
        PT_0B = [0, 0, 0]
        T_0B = [0, 0, 0]
        PT_AB = [0, 0, 0]
        T_AB = [0, 0, 0]
        PT_BB = [0, 0, 0]
        T_BB = [0, 0, 0]
        PT_ABB = [0, 0, 0]
        T_ABB = [0, 0, 0]
        F_BB = [0, 0, 0]
        F_ABB = [0, 0, 0]
        Succ_A = [0, 0, 0]
        Succ_B = [0, 0, 0]
        Succ_M = [0, 0, 0]
        Succ_back = [0, 0, 0]
        finalfail = [0, 0, 0]

        change = [0, 0, 0]
        for i in range(len(self.transfer)):
            ''''
            if self.stage == 0:
                np.random.seed(i)
            '''
            for j in range(3):

                if i == 0: #dynamic transfer value for new infection with different layers j
                    self.transfer[0][1] = self.s_0[j] * (1 - self.alpha_0)
                    self.transfer[0][2] = self.s_0[j] * self.alpha_0
                    self.transfer[0][3] = self.s_a[j] * (1 - self.alpha_0)
                    self.transfer[0][4] = self.s_a[j] * self.alpha_0
                    self.transfer[0][5] = self.s_b[j] * (1 - self.alpha_0)
                    self.transfer[0][6] = self.s_b[j] * self.alpha_0
                    self.transfer[0][7] = self.s_ab[j] * (1 - self.alpha_0)
                    self.transfer[0][8] = self.s_ab[j] * self.alpha_0
                    self.transfer[0][0] = 1 - sum(self.transfer[0][1:])
                #change = self.rng.multinomial(self.init_state[i], self.transfer[i])

                change[j] = self.rng.multinomial(self.init_state[j][i], self.transfer[i])
                #print('initial_state: ', self.init_state)
                #print('transfer: ', self.transfer)
                self.state[j] += change[j]


                '''
                if i == 0:
                    inci_in += sum(change[1:9])
                    N_0 += max(0, change[1])
                    I_0 += max(0, change[2])
                    N_A += max(0, change[3])
                    I_A += max(0, change[4])
                    N_B += max(0, change[5])
                    I_B += max(0, change[6])
                    N_AB += max(0, change[7])
                    I_AB += max(0, change[8])
                    QALYlost += self.UR * self.JR * (max(0, change[2])+max(0, change[4])+max(0, change[6])+max(0, change[8]))
                    Infection_qalyloss += self.UR * self.JR * (max(0, change[2])+max(0, change[4])+max(0, change[6])+max(0, change[8]))
    
                if i == 1:
                    PT0_in += max(0, change[9]) + max(0, change[19])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[9]) + (self.CB + self.Diagnose) * max(0, change[19]) + self.CS * max(0, change[0])
                    QALYlost += self.QS * max(0, change[0])
                    PT_0A += max(0, change[9])
                    PT_0B += max(0, change[19])
    
                    A_cost += self.CA * max(0, change[9])
                    B_cost += self.CB * max(0, change[19])
                    Diagnose_cost += self.Diagnose * max(0, change[9])+ self.Diagnose * max(0, change[19])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
    
                if i == 2:
                    T0_in += max(0, change[10]) + max(0, change[20])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[10]) + (self.CB + self.Diagnose) * max(0, change[20]) + self.CS * max(0, change[0])
                    QALYlost += self.UR * self.JA * max(0, change[10]) + self.UR * self.JB * max(0, change[20]) + self.QS * max(0, change[0])
                    T_0A += max(0, change[10])
                    T_0B += max(0, change[20])
    
                    A_cost += self.CA * max(0, change[10])
                    B_cost += self.CB * max(0, change[20])
                    Diagnose_cost += self.Diagnose * max(0, change[10])+ self.Diagnose * max(0, change[20])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    A_qalysloss += self.UR * self.JA * max(0, change[10])
                    B_qalysloss += self.UR * self.JB * max(0, change[20])
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
                if i == 3:
                    PTA_in += max(0, change[11]) + max(0, change[21])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[11]) + (self.CB + self.Diagnose) * max(0, change[21]) + self.CS * max(0, change[0])
                    QALYlost += self.QS * max(0, change[0])
                    PT_AA += max(0, change[11])
                    PT_AB += max(0, change[21])
    
                    A_cost += self.CA * max(0, change[11])
                    B_cost += self.CB * max(0, change[21])
                    Diagnose_cost += self.Diagnose * max(0, change[11])+ self.Diagnose * max(0, change[21])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
                if i == 4:
                    TA_in += max(0, change[12]) + max(0, change[22])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[12]) + (self.CB + self.Diagnose) * max(0, change[22]) + self.CS * max(0, change[0])
                    QALYlost += self.UR * self.JA * max(0, change[12]) + self.UR * self.JB * max(0, change[22]) + self.QS * max(0, change[0])
                    T_AA += max(0, change[12])
                    T_AB += max(0, change[22])
    
                    A_cost += self.CA * max(0, change[12])
                    B_cost += self.CB * max(0, change[22])
                    Diagnose_cost += self.Diagnose * max(0, change[12])+ self.Diagnose * max(0, change[22])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    A_qalysloss += self.UR * self.JA * max(0, change[12])
                    B_qalysloss += self.UR * self.JB * max(0, change[22])
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
    
                if i == 5:
                    PTB_in += max(0, change[13]) + max(0, change[23])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[13]) + (self.CB + self.Diagnose) * max(0, change[23]) + self.CS * max(0, change[0])
                    QALYlost += self.QS * max(0, change[0])
                    PT_BA += max(0, change[13])
                    PT_BB += max(0, change[23])
    
                    A_cost += self.CA * max(0, change[13])
                    B_cost += self.CB * max(0, change[23])
                    Diagnose_cost += self.Diagnose * max(0, change[13])+ self.Diagnose * max(0, change[23])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
                if i == 6:
                    TB_in += max(0, change[14]) + max(0, change[24])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[14]) + (self.CB + self.Diagnose) * max(0, change[24]) + self.CS * max(0, change[0])
                    QALYlost += self.UR * self.JA * max(0, change[14]) + self.UR * self.JB * max(0, change[24]) + self.QS * max(0, change[0])
                    T_BA += max(0, change[14])
                    T_BB += max(0, change[24])
    
                    A_cost += self.CA * max(0, change[14])
                    B_cost += self.CB * max(0, change[24])
                    Diagnose_cost += self.Diagnose * max(0, change[14])+ self.Diagnose * max(0, change[24])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    A_qalysloss += self.UR * self.JA * max(0, change[14])
                    B_qalysloss += self.UR * self.JB * max(0, change[24])
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
                if i == 7:
                    PTAB_in += max(0, change[15]) + max(0, change[25])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[15]) + (self.CB + self.Diagnose) * max(0, change[25]) + self.CS * max(0, change[0])
                    QALYlost += self.QS * max(0, change[0])
                    ABcost += self.CA * max(0, change[15]) + self.CB * max(0, change[25])
                    PT_ABA += max(0, change[15])
                    PT_ABB += max(0, change[25])
    
                    A_cost += self.CA * max(0, change[15])
                    B_cost += self.CB * max(0, change[25])
                    Diagnose_cost += self.Diagnose * max(0, change[15])+ self.Diagnose * max(0, change[25])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
                if i == 8:
                    TAB_in += max(0, change[16]) + max(0, change[26])
                    Totalcost += (self.CA + self.Diagnose) * max(0, change[16]) + (self.CB + self.Diagnose) * max(0, change[26]) + self.CS * max(0, change[0])
                    QALYlost += self.UR * self.JA * max(0, change[16]) + self.UR * self.JB * max(0, change[26]) + self.QS * max(0, change[0])
                    ABcost += (self.CA + self.Diagnose) * max(0, change[16]) + (self.CB + self.Diagnose) * max(0, change[26])
                    ABQALY += self.UR * self.JA * max(0, change[16]) + self.UR * self.JB * max(0, change[26])
                    AfterswitchAB += max(0, change[26])
                    T_ABA += max(0, change[16])
                    T_ABB += max(0, change[26])
    
                    A_cost += self.CA * max(0, change[16])
                    B_cost += self.CB * max(0, change[26])
                    Diagnose_cost += self.Diagnose * max(0, change[16])+ self.Diagnose * max(0, change[26])
                    Sequalae_cost += self.CS * max(0, change[0])
    
                    A_qalysloss += self.UR * self.JA * max(0, change[16])
                    B_qalysloss += self.UR * self.JB * max(0, change[26])
                    Sequalae_qalysloss += self.QS * max(0, change[0])
    
                if i == 9:
                    PT_AA += max(0, change[11])
                    Succ_A += max(0, change[0])
    
                if i == 10:
                    T_AA += max(0, change[12])
                    Succ_A += max(0, change[0])
    
                if i == 11:
                    N_A += max(0, change[3])
    
    
                if i == 12:
                    Totalcost += self.CS * max(0, change[17])
                    QALYlost += self.QS * max(0, change[17])
                    BswitchAB += max(0, change[17])
                    F_A += max(0, change[17])
    
                    Sequalae_cost += self.CS * max(0, change[17])
    
                    Sequalae_qalysloss += self.QS * max(0, change[17])
    
                if i == 13:
                    PT_ABA += max(0, change[15])
                    Succ_A += max(0, change[0])
    
                if i == 14:
                    T_ABA += max(0, change[16])
                    Succ_A += max(0, change[0])
    
                if i == 15:
                    N_AB += max(0, change[7])
    
                if i == 16:
                    Totalcost += self.CS * max(0, change[18])
                    QALYlost += self.QS * max(0, change[18])
                    ABcost += self.CS * max(0, change[18])
                    ABQALY += self.QS * max(0, change[18])
                    BswitchAB += max(0, change[18])
                    F_ABA += max(0, change[18])
    
                    Sequalae_cost += self.CS * max(0, change[18])
    
                    Sequalae_qalysloss += self.QS * max(0, change[18])
    
                if i == 17:
                    Totalcost += self.CB * max(0, change[22])
                    QALYlost += self.UR * self.JB * max(0, change[22])
                    FailAcost += self.CB * max(0, change[22])
                    FailAQALY += self.UR * self.JB * max(0, change[22])
                    BswitchAB += max(0, change[22])
                    T_AB += max(0, change[22])
    
                    B_cost += self.CB * max(0, change[22])
    
                    B_qalysloss += self.UR * self.JB * max(0, change[22])
    
                if i == 18:
                    Totalcost += self.CB * max(0, change[26])
                    QALYlost += self.UR * self.JB * max(0, change[26])
                    FailAcost += self.CB * max(0, change[26])
                    FailAQALY += self.UR * self.JB * max(0, change[26])
                    ABcost += self.CB * max(0, change[26])
                    ABQALY += self.UR * self.JB * max(0, change[26])
                    BeforeswitchAB += max(0, change[26])
                    AfterswitchAB += max(0, change[26])
                    BswitchAB += max(0, change[26])
                    T_ABB += max(0, change[26])
    
                    B_cost += self.CB * max(0, change[26])
    
                    B_qalysloss += self.UR * self.JB * max(0, change[26])
    
                if i == 19:
                    PT_BB += max(0, change[23])
                    Succ_B += max(0, change[0])
    
                if i == 20:
                    T_BB += max(0, change[24])
                    Succ_B += max(0, change[0])
    
                if i == 21:
                    PT_ABB += max(0, change[25])
                    Succ_B += max(0, change[0])
    
                if i == 22:
                    BeforeswitchAB += max(0, change[26])
                    AfterswitchAB += max(0, change[26])
                    T_ABB += max(0, change[26])
                    Succ_B += max(0, change[0])
    
                if i == 23:
                    N_B += max(0, change[5])
    
    
                if i == 24:
                    Totalcost += self.CS * max(0, change[27])
                    QALYlost += self.QS * max(0, change[27])
                    BswitchAB += max(0, change[27])
                    F_BB += max(0, change[27])
    
                    Sequalae_cost += self.CS * max(0, change[27])
    
                    Sequalae_qalysloss += self.QS * max(0, change[27])
    
                if i == 25:
                    N_AB += max(0, change[7])
    
                if i == 26:
                    Totalcost += self.CS * max(0, change[28])
                    QALYlost += self.QS * max(0, change[28])
                    ABcost += self.CS * max(0, change[28])
                    ABQALY += self.QS * max(0, change[28])
                    BswitchAB += max(0, change[28])
                    F_ABB += max(0, change[28])
    
                    Sequalae_cost += self.CS * max(0, change[28])
    
                    Sequalae_qalysloss += self.QS * max(0, change[28])
    
                if i == 27:
                    Totalcost += (self.CA + self.CD) * max(0, change[0])
                    QALYlost += self.UR * self.JA  * max(0, change[0])
                    FailBcost += (self.CA + self.CD) * max(0, change[0])
                    FailBQALY += self.UR * self.JA  * max(0, change[0])
                    BswitchAB += max(0, change[0])
                    Succ_A += max(0, change[0])
    
                    A_cost += self.CA * max(0, change[0])
    
                    Drugtest_cost += self.CD * max(0, change[0])
    
    
    
                if i == 28:
                    #Totalcost += (self.CM + self.CD) * max(0, change[0])
                    #QALYlost += self.UM * self.JM * max(0, change[0])
                    FailBcost += (self.CM + self.CD) * max(0, change[0])
                    FailBQALY += self.UM * self.JM * max(0, change[0])
                    ABcost += (self.CM + self.CD) * max(0, change[0])
                    ABQALY += self.UM * self.JM * max(0, change[0])
                    BswitchAB += max(0, change[0])
                    if not self.M_flag:
                        Succ_M += max(0, change[0])
                        Totalcost += (self.CM + self.CD) * max(0, change[0])
                        QALYlost += self.UM * self.JM * max(0, change[0])
    
                        M_cost += self.CM * max(0, change[0])
                        Drugtest_cost += self.CD * max(0, change[0])
    
                        M_qalysloss += self.UM * self.JM * max(0, change[0])
    
                    else:
                        Succ_back += max(0, change[0])
                        Totalcost += (self.CB + self.CD) * max(0, change[0])
                        QALYlost += self.UR * self.JB * max(0, change[0])
    
                        B_cost += self.CB * max(0, change[0])
                        Drugtest_cost += self.CD * max(0, change[0])
    
                        B_qalysloss += self.UR * self.JB * max(0, change[26])
    
                    TC_in += max(0, change[29])
    
    
                if i == 29:
                    Totalcost += self.CM * max(0, change[0])
                    QALYlost += self.UM * self.JM * max(0, change[0])
                    Succ_M += max(0, change[0])
    
                    M_cost += self.CM * max(0, change[0])
    
                    M_qalysloss += self.UM * self.JM * max(0, change[0])
                '''#old CEA
                if i == 0:
                    inci_in[j] += sum(change[j][1:9])
                    N_0[j] += max(0, change[j][1])
                    I_0[j] += max(0, change[j][2])
                    N_A[j] += max(0, change[j][3])
                    I_A[j] += max(0, change[j][4])
                    N_B[j] += max(0, change[j][5])
                    I_B[j] += max(0, change[j][6])
                    N_AB[j] += max(0, change[j][7])
                    I_AB[j] += max(0, change[j][8])
                    Totalcost[j] += self.CS * (max(0, change[j][1])+max(0, change[j][3])+max(0, change[j][5])+max(0, change[j][7]))
                    Sequalae_cost[j] += self.CS * (max(0, change[j][1])+max(0, change[j][3])+max(0, change[j][5])+max(0, change[j][7]))
                    QALYlost[j] += self.UR * self.JR * (max(0, change[j][2])+max(0, change[j][4])+max(0, change[j][6])+max(0, change[j][8])) + self.QS * (max(0, change[j][1])+max(0, change[j][3])+max(0, change[j][5])+max(0, change[j][7]))
                    Infection_qalyloss[j] += self.UR * self.JR * (max(0, change[j][2])+max(0, change[j][4])+max(0, change[j][6])+max(0, change[j][8]))
                    Sequalae_qalysloss[j] += self.QS * (max(0, change[j][1])+max(0, change[j][3])+max(0, change[j][5])+max(0, change[j][7]))

                if i == 1:
                    PT0_in[j] += max(0, change[j][9]) + max(0, change[j][19])
                    Totalcost[j] += (self.CA2 + self.Diagnose) * max(0, change[j][9]) + (self.CB2 + self.Diagnose) * max(0, change[j][19])
                    PT_0A[j] += max(0, change[j][9])
                    PT_0B[j] += max(0, change[j][19])

                    A_cost[j] += self.CA2 * max(0, change[j][9])
                    B_cost[j] += self.CB2 * max(0, change[j][19])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][9])+ self.Diagnose * max(0, change[j][19])



                if i == 2:
                    T0_in[j] += max(0, change[j][10]) + max(0, change[j][20])
                    Totalcost[j] += (self.CA1 + self.Diagnose) * max(0, change[j][10]) + (self.CB1 + self.Diagnose) * max(0, change[j][20]) + self.CS * max(0, change[j][0])
                    QALYlost[j] += self.QS * max(0, change[j][0])
                    T_0A[j] += max(0, change[j][10])
                    T_0B[j] += max(0, change[j][20])

                    A_cost[j] += self.CA1 * max(0, change[j][10])
                    B_cost[j] += self.CB1 * max(0, change[j][20])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][10])+ self.Diagnose * max(0, change[j][20])
                    Sequalae_cost[j] += self.CS * max(0, change[j][0])

                    Sequalae_qalysloss[j] += self.QS * max(0, change[j][0])

                if i == 3:
                    PTA_in[j] += max(0, change[j][11]) + max(0, change[j][21])
                    Totalcost[j] += (self.CA2 + self.Diagnose) * max(0, change[j][11]) + (self.CB2 + self.Diagnose) * max(0, change[j][21])
                    PT_AA[j] += max(0, change[j][11])
                    PT_AB[j] += max(0, change[j][21])

                    A_cost[j] += self.CA2 * max(0, change[j][11])
                    B_cost[j] += self.CB2 * max(0, change[j][21])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][11])+ self.Diagnose * max(0, change[j][21])


                if i == 4:
                    TA_in[j] += max(0, change[j][12]) + max(0, change[j][22])
                    Totalcost[j] += (self.CA1 + self.Diagnose) * max(0, change[j][12]) + (self.CB1 + self.Diagnose) * max(0, change[j][22]) + self.CS * (max(0, change[j][0])+max(0, change[j][12]))
                    QALYlost[j] += self.QS * (max(0, change[j][0])+max(0, change[j][12]))
                    T_AA[j] += max(0, change[j][12])
                    T_AB[j] += max(0, change[j][22])

                    A_cost[j] += self.CA1 * max(0, change[j][12])
                    B_cost[j] += self.CB1 * max(0, change[j][22])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][12])+ self.Diagnose * max(0, change[j][22])
                    Sequalae_cost[j] += self.CS * (max(0, change[j][0])+max(0, change[j][12]))

                    Sequalae_qalysloss[j] += self.QS * (max(0, change[j][0])+max(0, change[j][12]))


                if i == 5:
                    PTB_in[j] += max(0, change[j][13]) + max(0, change[j][23])
                    Totalcost[j] += (self.CA2 + self.Diagnose) * max(0, change[j][13]) + (self.CB2 + self.Diagnose) * max(0, change[j][23])
                    PT_BA[j] += max(0, change[j][13])
                    PT_BB[j] += max(0, change[j][23])

                    A_cost[j] += self.CA2 * max(0, change[j][13])
                    B_cost[j] += self.CB2 * max(0, change[j][23])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][13])+ self.Diagnose * max(0, change[j][23])


                if i == 6:
                    TB_in[j] += max(0, change[j][14]) + max(0, change[j][24])
                    Totalcost[j] += (self.CA1 + self.Diagnose) * max(0, change[j][14]) + (self.CB1 + self.Diagnose) * max(0, change[j][24]) + self.CS * (max(0, change[j][0]) + max(0, change[j][24]))
                    QALYlost[j] += self.QS * (max(0, change[j][0]) + max(0, change[j][24]))
                    T_BA[j] += max(0, change[j][14])
                    T_BB[j] += max(0, change[j][24])

                    A_cost[j] += self.CA1 * max(0, change[j][14])
                    B_cost[j] += self.CB1 * max(0, change[j][24])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][14])+ self.Diagnose * max(0, change[j][24])
                    Sequalae_cost[j] += self.CS * (max(0, change[j][0]) + max(0, change[j][24]))

                    Sequalae_qalysloss[j] += self.QS * (max(0, change[j][0]) + max(0, change[j][24]))

                if i == 7:
                    PTAB_in[j] += max(0, change[j][15]) + max(0, change[j][25])
                    Totalcost[j] += (self.CA2 + self.Diagnose) * max(0, change[j][15]) + (self.CB2 + self.Diagnose) * max(0, change[j][25])
                    ABcost[j] += self.CA2 * max(0, change[j][15]) + self.CB2 * max(0, change[j][25])
                    PT_ABA[j] += max(0, change[j][15])
                    PT_ABB[j] += max(0, change[j][25])

                    A_cost[j] += self.CA2 * max(0, change[j][15])
                    B_cost[j] += self.CB2 * max(0, change[j][25])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][15])+ self.Diagnose * max(0, change[j][25])

                if i == 8:
                    TAB_in[j] += max(0, change[j][16]) + max(0, change[j][26])
                    Totalcost[j] += (self.CA1 + self.Diagnose) * max(0, change[j][16]) + (self.CB1 + self.Diagnose) * max(0, change[j][26]) + self.CS * (max(0, change[j][0]) + max(0, change[j][16]) + max(0, change[j][26]))
                    QALYlost[j] += self.QS * (max(0, change[j][0]) + max(0, change[j][16]) + max(0, change[j][26]))
                    ABcost[j] += (self.CA1 + self.Diagnose) * max(0, change[j][16]) + (self.CB1 + self.Diagnose) * max(0, change[j][26])
                    AfterswitchAB[j] += max(0, change[j][26])
                    T_ABA[j] += max(0, change[j][16])
                    T_ABB[j] += max(0, change[j][26])

                    A_cost[j] += self.CA1 * max(0, change[j][16])
                    B_cost[j] += self.CB1 * max(0, change[j][26])
                    Diagnose_cost[j] += self.Diagnose * max(0, change[j][16])+ self.Diagnose * max(0, change[j][26])
                    Sequalae_cost[j] += self.CS * (max(0, change[j][0]) + max(0, change[j][16]) + max(0, change[j][26]))

                    Sequalae_qalysloss[j] += self.QS * (max(0, change[j][0]) + max(0, change[j][16]) + max(0, change[j][26]))

                if i == 9:
                    PT_AA[j] += max(0, change[j][11])
                    Succ_A[j] += max(0, change[j][0])

                if i == 10:
                    T_AA[j] += max(0, change[j][12])
                    Succ_A[j] += max(0, change[j][0])

                if i == 11:
                    N_A[j] += max(0, change[j][3])


                if i == 12:
                    BswitchAB[j] += max(0, change[j][17])
                    F_A[j] += max(0, change[j][17])
                    treatment_fail[j] += max(0, change[j][17])


                if i == 13:
                    PT_ABA[j] += max(0, change[j][15])
                    Succ_A[j] += max(0, change[j][0])

                if i == 14:
                    T_ABA[j] += max(0, change[j][16])
                    Succ_A[j] += max(0, change[j][0])

                if i == 15:
                    N_AB[j] += max(0, change[j][7])

                if i == 16:
                    BswitchAB[j] += max(0, change[j][18])
                    F_ABA[j] += max(0, change[j][18])
                    treatment_fail[j] += max(0, change[j][18])


                if i == 17:
                    Totalcost[j] += self.CB2 * max(0, change[j][22])
                    FailAcost[j] += self.CB2 * max(0, change[j][22])
                    BswitchAB[j] += max(0, change[j][22])
                    T_AB[j] += max(0, change[j][22])

                    B_cost[j] += self.CB2 * max(0, change[j][22])


                if i == 18:
                    Totalcost[j] += self.CB2 * max(0, change[j][26])
                    FailAcost[j] += self.CB2 * max(0, change[j][26])
                    ABcost[j] += self.CB2 * max(0, change[j][26])
                    BeforeswitchAB[j] += max(0, change[j][26])
                    AfterswitchAB[j] += max(0, change[j][26])
                    BswitchAB[j] += max(0, change[j][26])
                    T_ABB[j] += max(0, change[j][26])

                    B_cost[j] += self.CB2 * max(0, change[j][26])


                if i == 19:
                    PT_BB[j] += max(0, change[j][23])
                    Succ_B[j] += max(0, change[j][0])

                if i == 20:
                    T_BB[j] += max(0, change[j][24])
                    Succ_B[j] += max(0, change[j][0])

                if i == 21:
                    PT_ABB[j] += max(0, change[j][25])
                    T_ABB[j] += max(0, change[j][26])
                    Succ_B[j] += max(0, change[j][0])

                if i == 22:
                    BeforeswitchAB[j] += max(0, change[j][26])
                    AfterswitchAB[j] += max(0, change[j][26])
                    T_ABB[j] += max(0, change[j][26])
                    Succ_B[j] += max(0, change[j][0])

                if i == 23:
                    N_B[j] += max(0, change[j][5])


                if i == 24:
                    BswitchAB[j] += max(0, change[j][27])
                    F_BB[j] += max(0, change[j][27])
                    treatment_fail[j] += max(0, change[j][27])


                if i == 25:
                    N_AB[j] += max(0, change[j][7])

                if i == 26:
                    BswitchAB[j] += max(0, change[j][28])
                    F_ABB[j] += max(0, change[j][28])
                    treatment_fail[j] += max(0, change[j][28])

                if i == 27:

                    if self.strategy_flag == 1: #with DSTs
                        Totalcost[j] += (self.CA2 + self.CD) * (max(0, change[j][0]) + max(0, change[j][28]))
                        FailBcost[j] += (self.CA2 + self.CD) * (max(0, change[j][0]) + max(0, change[j][28]))
                        #F_ABB += max(0, change[j][28])
                        BswitchAB[j] += max(0, change[j][0])
                        Succ_A[j] += max(0, change[j][0]) + max(0, change[j][28])

                        finalfail[j] += max(0, change[j][28])

                        A_cost[j] += self.CA2 * max(0, change[j][0])

                        Drugtest_cost[j] += self.CD * max(0, change[j][0])

                        DST_in[j] += max(0, change[j][0])

                    else: #no DSTs
                        if not self.M_flag:
                            Totalcost[j] += self.CM * max(0, change[j][0])
                            QALYlost[j] += self.UM * self.JM * max(0, change[j][0])
                            FailBcost[j] += self.CM * max(0, change[j][0])
                            BswitchAB[j] += max(0, change[j][0])

                            Succ_M[j] += max(0, change[j][0])
                            M_cost[j] += self.CM * max(0, change[j][0])
                            M_qalysloss[j] += self.UM * self.JM * max(0, change[j][0])
                        else:
                            FailBcost[j] += self.CB2 * max(0, change[j][0])
                            ABcost[j] += self.CB2 * max(0, change[j][0])
                            Succ_back[j] += max(0, change[j][0])
                            Totalcost[j] += self.CB2 * max(0, change[j][0])

                            B_cost[j] += self.CB2 * max(0, change[j][0])
                            Succ_B[j] += max(0, change[j][0])

                            TC_in[j] += max(0, change[j][29])
                            treatment_fail[j] += max(0, change[j][29])


                if i == 28:
                    if self.strategy_flag == 1: #DSTs
                        BswitchAB[j] += max(0, change[j][0])
                        if not self.M_flag:
                            FailBcost[j] += (self.CM + self.CD) * max(0, change[j][0])
                            FailBQALY[j] += self.UM * self.JM * max(0, change[j][0])
                            ABcost[j] += (self.CM + self.CD) * max(0, change[j][0])
                            ABQALY[j] += self.UM * self.JM * max(0, change[j][0])
                            Succ_M[j] += max(0, change[j][0])
                            Totalcost[j] += (self.CM + self.CD) * max(0, change[j][0])
                            QALYlost[j] += self.UM * self.JM * max(0, change[j][0])

                            M_cost[j] += self.CM * max(0, change[j][0])
                            Drugtest_cost[j] += self.CD * (max(0, change[j][0]) + max(0, change[j][29]))

                            M_qalysloss[j] += self.UM * self.JM * max(0, change[j][0])

                        else:
                            FailBcost[j] += (self.CB2 + self.CD) * max(0, change[j][0])
                            ABcost[j] += (self.CB2 + self.CD) * max(0, change[j][0])
                            Succ_back[j] += max(0, change[j][0])
                            Totalcost[j] += (self.CB2 + self.CD) * max(0, change[j][0])

                            B_cost[j] += self.CB2 * max(0, change[j][0])
                            Drugtest_cost[j] += self.CD * (max(0, change[j][0]) + max(0, change[j][29]))


                        TC_in[j] += max(0, change[j][29])
                        treatment_fail[j] += max(0, change[j][29])
                        DST_in[j] += max(0, change[j][29]) + max(0, change[j][0]) #remove when using M

                    else:  #no DSTs
                        BswitchAB[j] += max(0, change[j][0])
                        if not self.M_flag:
                            FailBcost[j] += self.CM * max(0, change[j][0])
                            FailBQALY[j] += self.UM * self.JM * max(0, change[j][0])
                            ABcost[j] += self.CM * max(0, change[j][0])
                            ABQALY[j] += self.UM * self.JM * max(0, change[j][0])
                            Succ_M[j] += max(0, change[j][0])
                            Totalcost[j] += self.CM * max(0, change[j][0])
                            QALYlost[j] += self.UM * self.JM * max(0, change[j][0])

                            M_cost[j] += self.CM * max(0, change[j][0])


                            M_qalysloss[j] += self.UM * self.JM * max(0, change[j][0])

                        else:
                            FailBcost[j] += self.CB2 * max(0, change[j][0])
                            ABcost[j] += self.CB2 * max(0, change[j][0])
                            Succ_back[j] += max(0, change[j][0])
                            Totalcost[j] += self.CB2 * max(0, change[j][0])

                            B_cost[j] += self.CB2 * max(0, change[j][0])
                            Succ_B[j] += max(0, change[j][0])


                            TC_in[j] += max(0, change[j][29])
                            treatment_fail[j] += max(0, change[j][29])




                if i == 29:
                    Totalcost[j] += self.CM * max(0, change[j][0])
                    QALYlost[j] += self.UM * self.JM * max(0, change[j][0])
                    Succ_M[j] += max(0, change[j][0])

                    M_cost[j] += self.CM * max(0, change[j][0])

                    M_qalysloss[j] += self.UM * self.JM * max(0, change[j][0])

            self.totalstate = sum(self.state[j] for j in range(3))

        return PT0_in, T0_in, PTA_in, TA_in, PTB_in, TB_in, PTAB_in, TAB_in, TC_in, inci_in, Totalcost, QALYlost, FailAcost, FailAQALY, FailBcost, FailBQALY, ABcost, ABQALY, BeforeswitchAB, AfterswitchAB, BswitchAB,\
               N_0, I_0, N_A, I_A, N_B, I_B, N_AB, I_AB, PT_0A, T_0A, PT_AA, T_AA, PT_BA, T_BA, PT_ABA, T_ABA, F_A, F_ABA, PT_0B, T_0B, PT_AB, T_AB, PT_BB, T_BB, PT_ABB, T_ABB, F_BB, F_ABB, Succ_A, Succ_B, Succ_M, Succ_back,\
               A_cost, B_cost, M_cost, Diagnose_cost, Drugtest_cost, Sequalae_cost, A_qalysloss, B_qalysloss, M_qalysloss, Sequalae_qalysloss, Infection_qalyloss, treatment_fail, DST_in, finalfail



    def print_matrix(self, matrix):
        for i in matrix:
            print('\t'.join([str(n) for n in i]))
        print('\n')

    def run(self):
        self.judge = 0
        cumu_PT0 = 0
        cumu_T0 = 0
        cumu_PTA = 0
        cumu_TA = 0
        cumu_PTB = 0
        cumu_TB = 0
        cumu_PTAB = 0
        cumu_TAB = 0
        cumu_TC = 0
        cumu_total = 0

        cumu_treatment_fail = 0
        cumu_DST_in = 0


        cumu_inci = 0
        self.stage = 0
        cumu_Asy = 0
        cumu_Sy = 0
        cumu_Totalcost = 0
        cumu_QALYlost = 0
        cumu_FAcost = 0
        cumu_FAQALY = 0
        cumu_FBcost = 0
        cumu_FBQALY = 0
        cumu_ABcost = 0
        cumu_ABQALY = 0
        cumu_BeforeswitchAB = 0
        cumu_AfterswitchAB = 0
        cumu_BswitchAB = 0


        cumu_A_cost = 0
        cumu_B_cost = 0
        cumu_M_cost = 0
        cumu_Diagnose_cost = 0
        cumu_Drugtest_cost = 0
        cumu_Sequalae_cost = 0

        cumu_A_qalysloss = 0
        cumu_B_qalysloss = 0
        cumu_M_qalysloss = 0
        cumu_Sequalae_qalysloss = 0
        cumu_Infection_qalyloss = 0




        cumu_N_0 = 0
        cumu_I_0 = 0
        cumu_N_A = 0
        cumu_I_A = 0
        cumu_N_B = 0
        cumu_I_B = 0
        cumu_N_AB = 0
        cumu_I_AB = 0
        cumu_PT_0A = 0
        cumu_T_0A = 0
        cumu_PT_AA = 0
        cumu_T_AA = 0
        cumu_PT_BA = 0
        cumu_T_BA = 0
        cumu_PT_ABA = 0
        cumu_T_ABA = 0
        cumu_F_A = 0
        cumu_F_ABA = 0
        cumu_PT_0B = 0
        cumu_T_0B = 0
        cumu_PT_AB = 0
        cumu_T_AB = 0
        cumu_PT_BB = 0
        cumu_T_BB = 0
        cumu_PT_ABB = 0
        cumu_T_ABB = 0
        cumu_F_BB = 0
        cumu_F_ABB = 0
        cumu_finalfail = 0
        cumu_Succ_A = 0
        cumu_Succ_B = 0
        cumu_Succ_M = 0
        cumu_Succ_back = 0

        cumu_annualCost = 0
        cumu_annualQALY = 0


        for week in range(2600):
            self.judge += 1
            ar = np.exp(-self.b_i * (week - self.t_0))
            #self.gamma = self.b_min + (1 - self.b_min) / (1 + ar)
            #print('week: ', week)
            #print('gamma: ', self.gamma)
            
            if week == 0:
                self.update_transfer()
            
            
            PT0_in, T0_in, PTA_in, TA_in, PTB_in, TB_in, PTAB_in, TAB_in, TC_in, inci_in, Totalcost, QALYlost, FailAcost, FailAQALY, FailBcost, FailBQALY, ABcost, ABQALY, BeforeswitchAB, AfterswitchAB, BswitchAB,\
                N_0, I_0, N_A, I_A, N_B, I_B, N_AB, I_AB, PT_0A, T_0A, PT_AA, T_AA, PT_BA, T_BA, PT_ABA, T_ABA, F_A, F_ABA, PT_0B, T_0B, PT_AB, T_AB, PT_BB, T_BB, PT_ABB, T_ABB, F_BB, F_ABB, Succ_A, Succ_B, Succ_M, Succ_back,\
                A_cost, B_cost, M_cost, Diagnose_cost, Drugtest_cost, Sequalae_cost, A_qalysloss, B_qalysloss, M_qalysloss, Sequalae_qalysloss, Infection_qalyloss, treatment_fail, DST_in, finalfail = self.update_state()




            cumu_PT0 += sum(PT0_in)
            cumu_T0 += sum(T0_in)
            cumu_PTA += sum(PTA_in)
            cumu_TA += sum(TA_in)
            cumu_PTB += sum(PTB_in)
            cumu_TB += sum(TB_in)
            cumu_PTAB += sum(PTAB_in)
            cumu_TAB += sum(TAB_in)
            cumu_TC += sum(TC_in)

            cumu_treatment_fail += sum(treatment_fail)
            cumu_DST_in += sum(DST_in)



            cumu_inci += sum(inci_in)

            cumu_Totalcost += sum(Totalcost)
            cumu_QALYlost += sum(QALYlost)


            cumu_A_cost += sum(A_cost)
            cumu_B_cost += sum(B_cost)
            cumu_M_cost += sum(M_cost)
            cumu_Diagnose_cost += sum(Diagnose_cost)
            cumu_Drugtest_cost += sum(Drugtest_cost)
            cumu_Sequalae_cost += sum(Sequalae_cost)

            cumu_A_qalysloss += sum(A_qalysloss)
            cumu_B_qalysloss += sum(B_qalysloss)
            cumu_M_qalysloss += sum(M_qalysloss)
            cumu_Sequalae_qalysloss += sum(Sequalae_qalysloss)
            cumu_Infection_qalyloss += sum(Infection_qalyloss)




            cumu_annualCost += sum(Totalcost)
            cumu_annualQALY += sum(QALYlost)

            cumu_FAcost += sum(FailAcost)
            cumu_FAQALY += sum(FailAQALY)
            cumu_FBcost += sum(FailBcost)
            cumu_FBQALY += sum(FailBQALY)
            cumu_ABcost += sum(ABcost)
            cumu_ABQALY += sum(ABQALY)
            cumu_BeforeswitchAB += sum(BeforeswitchAB)
            cumu_AfterswitchAB += sum(AfterswitchAB)
            cumu_BswitchAB += sum(BswitchAB)


            cumu_N_0 += sum(N_0)
            cumu_I_0 += sum(I_0)
            cumu_N_A += sum(N_A)
            cumu_I_A += sum(I_A)
            cumu_N_B += sum(N_B)
            cumu_I_B += sum(I_B)
            cumu_N_AB += sum(N_AB)
            cumu_I_AB += sum(I_AB)
            cumu_PT_0A += sum(PT_0A)
            cumu_T_0A += sum(T_0A)
            cumu_PT_AA += sum(PT_AA)
            cumu_T_AA += sum(T_AA)
            cumu_PT_BA += sum(PT_BA)
            cumu_T_BA += sum(T_BA)
            cumu_PT_ABA += sum(PT_ABA)
            cumu_T_ABA += sum(T_ABA)
            cumu_F_A += sum(F_A)
            cumu_F_ABA += sum(F_ABA)
            cumu_PT_0B += sum(PT_0B)
            cumu_T_0B += sum(T_0B)
            cumu_PT_AB += sum(PT_AB)
            cumu_T_AB += sum(T_AB)
            cumu_PT_BB += sum(PT_BB)
            cumu_T_BB += sum(T_BB)
            cumu_PT_ABB += sum(PT_ABB)
            cumu_T_ABB += sum(T_ABB)
            cumu_F_BB += sum(F_BB)
            cumu_F_ABB += sum(F_ABB)
            cumu_Succ_A += sum(Succ_A)
            cumu_Succ_B += sum(Succ_B)
            cumu_Succ_M += sum(Succ_M)
            cumu_Succ_back += sum(Succ_back)
            cumu_finalfail += sum(finalfail)



            '''
            if stage == 0:
                cumu_T0 += T0_in
                cumu_TA += TA_in
                cumu_inci += inci_in
                cumu_Asy += Asy_in
                cumu_Sy += Sy_in
            else:
                cumu_TA = 0
            '''
            '''
            if self.state[2] + self.state[4] == 0:
                # failed simu
                return False
            '''
            #cur_resistant_rate = (self.state[2] + self.state[4]) / sum(self.state)
            Acur_resistant_rate = (self.totalstate[3] + self.totalstate[4]) / sum(self.totalstate[1:9])
            Bcur_resistant_rate = (self.totalstate[5] + self.totalstate[6]) / sum(self.totalstate[1:9])
            ABcur_resistant_rate = (self.totalstate[7] + self.totalstate[8]) / sum(self.totalstate[1:9])

            cur_prevalence = sum(self.totalstate[1:9]) / sum(self.totalstate)
            cur_symptomatic = (self.totalstate[2] + self.totalstate[4] + self.totalstate[6] + self.totalstate[8]) / sum(self.totalstate[1:9])
            #print('cur_prevalence: ', cur_prevalence)
            self.resistant_rate.append(Acur_resistant_rate)
            self.prevalence.append(cur_prevalence)
            self.symptomatic.append(cur_symptomatic)
            self.S.append(self.totalstate[0])
            self.AI.append(self.totalstate[1])
            self.AR.append(self.totalstate[2])
            self.I_0.append(self.totalstate[3])
            self.R_0.append(self.totalstate[4])
            self.T_I.append(self.totalstate[5])
            self.T_A.append(self.totalstate[6])
            self.T_M.append(self.totalstate[7])
            self.F_A.append(self.totalstate[8])

            '''
            if cur_prevalence < 0.01 or cur_prevalence > 0.15:
                # failed simu
                return False
            '''
            #print('cumu_T0: ', cumu_T0)
            #print('cumu_TA: ', cumu_TA)
            '''
            print(self.state, sum(self.state))
            print(f'Resistant Proportion: {(self.state[2] + self.state[4]) / sum(self.state[1:5])}')
            print(f'Prevalence: {(self.state[1] + self.state[2] + self.state[3] + self.state[4]) / sum(self.state)}')
            print(f"Cumu T0 and TA: {cumu_T0, cumu_TA}")
            '''
            #False and
            '''
            if self.stage ==0 and week == 2599 and (cumu_TA + cumu_PTA) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB) < 0.05:
                # failed simu
                return False


            if week == 2599 and (cumu_TB + cumu_PTB) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB) < 0.05:
                # failed simu
                return False
            '''


            if week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 0 and \
                    (cumu_TA + cumu_PTA) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.transfer[1][19] += self.transfer[1][9]
                self.transfer[2][20] += self.transfer[2][10]
                self.transfer[3][21] += self.transfer[3][11]
                self.transfer[4][22] += self.transfer[4][12]
                self.transfer[5][23] += self.transfer[5][13]
                self.transfer[6][24] += self.transfer[6][14]
                self.transfer[7][25] += self.transfer[7][15]
                self.transfer[8][26] += self.transfer[8][16]

                self.transfer[1][9] = 0
                self.transfer[2][10] = 0
                self.transfer[3][11] = 0
                self.transfer[4][12] = 0
                self.transfer[5][13] = 0
                self.transfer[6][14] = 0
                self.transfer[7][15] = 0
                self.transfer[8][16] = 0

                self.stage = 1

            #new drug for every n years



            if week == self.newdrug_time:
                self.M_flag = True
                self.transfer[28][29] = self.transfer[20][24]
                if self.strategy_flag == 0:
                    self.transfer[27][29] = self.transfer[20][24] #no drug test, use M for all last treatment case
                #self.CM = self.CA
                #self.UM = self.UR



            if week == self.newdrug_time_2:
                self.M_flag = True
                self.transfer[28][29] = self.transfer[20][24]
                if self.strategy_flag == 0:
                    self.transfer[27][29] = self.transfer[20][24]
                #self.CM = self.CA
                #self.UM = self.UR



            if week == self.newdrug_time_3:
                self.M_flag = True
                self.transfer[28][29] = self.transfer[20][24]
                if self.strategy_flag == 0:
                    self.transfer[27][29] = self.transfer[20][24]
                #self.CM = self.CA
                #self.UM = self.UR


            '''
            if week == 2559:
                self.M_flag = True
                self.transfer[28][29] = self.transfer[10][12]
                self.CM = self.CA
                self.UM = self.UR
            '''


             #ResistanceB (self.comeyear * 52 - 1)

            if week >= self.newdrug_time and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 1 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                for j in range(3):
                    self.state[j][1] += self.state[j][3]
                    self.state[j][2] += self.state[j][4]
                    self.state[j][3] = self.state[j][5] + self.state[j][7]
                    self.state[j][4] = self.state[j][6] + self.state[j][8]
                    self.state[j][5] = 0
                    self.state[j][6] = 0
                    self.state[j][7] = 0
                    self.state[j][8] = 0

                    self.state[j][9] += self.state[j][19] + self.state[j][21] + self.state[j][11]
                    self.state[j][10] += self.state[j][20] + self.state[j][22] + self.state[j][17] + self.state[j][12]
                    self.state[j][11] = self.state[j][23] + self.state[j][25] + self.state[j][15]
                    self.state[j][12] = self.state[j][24] + self.state[j][26] + self.state[j][18] + self.state[j][16]
                    #self.state[13] = 0
                    #self.state[14] = 0
                    self.state[j][15] = 0
                    self.state[j][16] = 0
                    self.state[j][17] = self.state[j][27] + self.state[j][28]
                    self.state[j][18] = 0

                    self.state[j][19] = 0
                    self.state[j][20] = 0
                    self.state[j][21] = 0
                    self.state[j][22] = 0
                    self.state[j][23] = 0
                    self.state[j][24] = 0
                    self.state[j][25] = 0
                    self.state[j][26] = 0

                    self.state[j][27] = 0
                    self.state[j][28] = 0

                if week < self.newdrug_time_2:
                    self.transfer[28][29] = 0
                    if self.strategy_flag == 0:
                        self.transfer[27][29] = 0 #no drug test, use M for all last treatment case
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 2
                self.M_flag = False


            if week >= self.newdrug_time_2 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 2 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                for j in range(3):
                    self.state[j][1] += self.state[j][3]
                    self.state[j][2] += self.state[j][4]
                    self.state[j][3] = self.state[j][5] + self.state[j][7]
                    self.state[j][4] = self.state[j][6] + self.state[j][8]
                    self.state[j][5] = 0
                    self.state[j][6] = 0
                    self.state[j][7] = 0
                    self.state[j][8] = 0

                    self.state[j][9] += self.state[j][19] + self.state[j][21] + self.state[j][11]
                    self.state[j][10] += self.state[j][20] + self.state[j][22] + self.state[j][17] + self.state[j][12]
                    self.state[j][11] = self.state[j][23] + self.state[j][25] + self.state[j][15]
                    self.state[j][12] = self.state[j][24] + self.state[j][26] + self.state[j][18] + self.state[j][16]
                    #self.state[j][13] = 0
                    #self.state[j][14] = 0
                    self.state[j][15] = 0
                    self.state[j][16] = 0
                    self.state[j][17] = self.state[j][27] + self.state[j][28]
                    self.state[j][18] = 0

                    self.state[j][19] = 0
                    self.state[j][20] = 0
                    self.state[j][21] = 0
                    self.state[j][22] = 0
                    self.state[j][23] = 0
                    self.state[j][24] = 0
                    self.state[j][25] = 0
                    self.state[j][26] = 0

                    self.state[j][27] = 0
                    self.state[j][28] = 0


                if week < self.newdrug_time_3:
                    self.transfer[28][29] = 0
                    if self.strategy_flag == 0:
                        self.transfer[27][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 3
                self.M_flag = False

            if week >= self.newdrug_time_3 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0)and self.stage == 3 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                for j in range(3):
                    self.state[j][1] += self.state[j][3]
                    self.state[j][2] += self.state[j][4]
                    self.state[j][3] = self.state[j][5] + self.state[j][7]
                    self.state[j][4] = self.state[j][6] + self.state[j][8]
                    self.state[j][5] = 0
                    self.state[j][6] = 0
                    self.state[j][7] = 0
                    self.state[j][8] = 0

                    self.state[j][9] += self.state[j][19] + self.state[j][21] + self.state[j][11]
                    self.state[j][10] += self.state[j][20] + self.state[j][22] + self.state[j][17] + self.state[j][12]
                    self.state[j][11] = self.state[j][23] + self.state[j][25] + self.state[j][15]
                    self.state[j][12] = self.state[j][24] + self.state[j][26] + self.state[j][18] + self.state[j][16]
                    #self.state[j][13] = 0
                    #self.state[j][14] = 0
                    self.state[j][15] = 0
                    self.state[j][16] = 0
                    self.state[j][17] = self.state[j][27] + self.state[j][28]
                    self.state[j][18] = 0

                    self.state[j][19] = 0
                    self.state[j][20] = 0
                    self.state[j][21] = 0
                    self.state[j][22] = 0
                    self.state[j][23] = 0
                    self.state[j][24] = 0
                    self.state[j][25] = 0
                    self.state[j][26] = 0

                    self.state[j][27] = 0
                    self.state[j][28] = 0


                if week < 2600:
                    self.transfer[28][29] = 0
                    if self.strategy_flag == 0:
                        self.transfer[27][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 4
                self.M_flag = False


            '''
            if week >= 519 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 1 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:

                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0
            
                if week < 1039:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 2
                self.M_flag = False


            
            if week >= 1039 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 2 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0



                if week < 1559:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 3
                self.M_flag = False


            if week >= 1559 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0)and self.stage == 3 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0


                if week < 2079:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 4
                self.M_flag = False
            

            if week >= 2079 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 4 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0


                if week < 2559:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 5
                self.M_flag = False


            if week >= 2599 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 5 and self.judge > 0 and \
                    (cumu_TB + cumu_PTB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0

                self.transfer[28][29] = 0
                self.CM = self.MoriCost
                self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 6
                self.M_flag = False
            '''



            #ResistanceAB

            if week >= self.newdrug_time and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 1 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                for j in range(3):
                    self.state[j][1] += self.state[j][3]
                    self.state[j][2] += self.state[j][4]
                    self.state[j][3] = self.state[j][5] + self.state[j][7]
                    self.state[j][4] = self.state[j][6] + self.state[j][8]
                    self.state[j][5] = 0
                    self.state[j][6] = 0
                    self.state[j][7] = 0
                    self.state[j][8] = 0

                    self.state[j][9] += self.state[j][19] + self.state[j][21] + self.state[j][11]
                    self.state[j][10] += self.state[j][20] + self.state[j][22] + self.state[j][17] + self.state[j][12]
                    self.state[j][11] = self.state[j][23] + self.state[j][25] + self.state[j][15]
                    self.state[j][12] = self.state[j][24] + self.state[j][26] + self.state[j][18] + self.state[j][16]
                    #self.state[j][13] = 0
                    #self.state[j][14] = 0
                    self.state[j][15] = 0
                    self.state[j][16] = 0
                    self.state[j][17] = self.state[j][27] + self.state[j][28]
                    self.state[j][18] = 0

                    self.state[j][19] = 0
                    self.state[j][20] = 0
                    self.state[j][21] = 0
                    self.state[j][22] = 0
                    self.state[j][23] = 0
                    self.state[j][24] = 0
                    self.state[j][25] = 0
                    self.state[j][26] = 0

                    self.state[j][27] = 0
                    self.state[j][28] = 0


                if week < self.newdrug_time_2:
                    self.transfer[28][29] = 0
                    if self.strategy_flag == 0:
                        self.transfer[27][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 2
                self.M_flag = False


            if week >= self.newdrug_time_2 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 2 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                for j in range(3):
                    self.state[j][1] += self.state[j][3]
                    self.state[j][2] += self.state[j][4]
                    self.state[j][3] = self.state[j][5] + self.state[j][7]
                    self.state[j][4] = self.state[j][6] + self.state[j][8]
                    self.state[j][5] = 0
                    self.state[j][6] = 0
                    self.state[j][7] = 0
                    self.state[j][8] = 0

                    self.state[j][9] += self.state[j][19] + self.state[j][21] + self.state[j][11]
                    self.state[j][10] += self.state[j][20] + self.state[j][22] + self.state[j][17] + self.state[j][12]
                    self.state[j][11] = self.state[j][23] + self.state[j][25] + self.state[j][15]
                    self.state[j][12] = self.state[j][24] + self.state[j][26] + self.state[j][18] + self.state[j][16]
                    #self.state[j][13] = 0
                    #self.state[j][14] = 0
                    self.state[j][15] = 0
                    self.state[j][16] = 0
                    self.state[j][17] = self.state[j][27] + self.state[j][28]
                    self.state[j][18] = 0

                    self.state[j][19] = 0
                    self.state[j][20] = 0
                    self.state[j][21] = 0
                    self.state[j][22] = 0
                    self.state[j][23] = 0
                    self.state[j][24] = 0
                    self.state[j][25] = 0
                    self.state[j][26] = 0

                    self.state[j][27] = 0
                    self.state[j][28] = 0



                if week < self.newdrug_time_3:
                    self.transfer[28][29] = 0
                    if self.strategy_flag == 0:
                        self.transfer[27][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 3
                self.M_flag = False


            if week >= self.newdrug_time_3 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0)and self.stage == 3 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                for j in range(3):
                    self.state[j][1] += self.state[j][3]
                    self.state[j][2] += self.state[j][4]
                    self.state[j][3] = self.state[j][5] + self.state[j][7]
                    self.state[j][4] = self.state[j][6] + self.state[j][8]
                    self.state[j][5] = 0
                    self.state[j][6] = 0
                    self.state[j][7] = 0
                    self.state[j][8] = 0

                    self.state[j][9] += self.state[j][19] + self.state[j][21] + self.state[j][11]
                    self.state[j][10] += self.state[j][20] + self.state[j][22] + self.state[j][17] + self.state[j][12]
                    self.state[j][11] = self.state[j][23] + self.state[j][25] + self.state[j][15]
                    self.state[j][12] = self.state[j][24] + self.state[j][26] + self.state[j][18] + self.state[j][16]
                    #self.state[j][13] = 0
                    #self.state[j][14] = 0
                    self.state[j][15] = 0
                    self.state[j][16] = 0
                    self.state[j][17] = self.state[j][27] + self.state[j][28]
                    self.state[j][18] = 0

                    self.state[j][19] = 0
                    self.state[j][20] = 0
                    self.state[j][21] = 0
                    self.state[j][22] = 0
                    self.state[j][23] = 0
                    self.state[j][24] = 0
                    self.state[j][25] = 0
                    self.state[j][26] = 0

                    self.state[j][27] = 0
                    self.state[j][28] = 0


                if week <= 2600:
                    self.transfer[28][29] = 0
                    if self.strategy_flag == 0:
                        self.transfer[27][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 4
                self.M_flag = False


            '''
            if week >= 519 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 1 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0
                
                
                if week <= 1039:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 2
                self.M_flag = False



            if week >= 1039 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 2 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0
                



                if week <= 1559:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 3
                self.M_flag = False


            if week >= 1559 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0)and self.stage == 3 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0


                if week <= 2079:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 4
                self.M_flag = False


            if week >= 2079 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 4 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0


                if week <= 2559:
                    self.transfer[28][29] = 0
                    self.CM = self.MoriCost
                    self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 5
                self.M_flag = False


            if week >= 2599 and week % 52 == 51 and (cumu_T0 + cumu_TA > 0) and self.stage == 5 and self.judge > 0 and \
                    (cumu_TAB + cumu_PTAB) / (cumu_T0 + cumu_TA + cumu_PT0 + cumu_PTA + cumu_TB + cumu_PTB + cumu_TAB + cumu_PTAB) > self.resistant_threshold:
                self.state[1] += self.state[3]
                self.state[2] += self.state[4]
                self.state[3] = self.state[5] + self.state[7]
                self.state[4] = self.state[6] + self.state[8]
                self.state[5] = 0
                self.state[6] = 0
                self.state[7] = 0
                self.state[8] = 0

                self.state[9] = self.state[19] + self.state[21]
                self.state[10] = self.state[20] + self.state[22]
                self.state[11] = self.state[23] + self.state[25]
                self.state[12] = self.state[24] + self.state[26]
                self.state[13] = 0
                self.state[14] = 0
                self.state[15] = 0
                self.state[16] = 0
                self.state[17] = self.state[27] + self.state[28]
                self.state[18] = 0

                self.state[19] = 0
                self.state[20] = 0
                self.state[21] = 0
                self.state[22] = 0
                self.state[23] = 0
                self.state[24] = 0
                self.state[25] = 0
                self.state[26] = 0

                self.state[27] = 0
                self.state[28] = 0

                self.transfer[28][29] = 0
                self.CM = self.MoriCost
                self.UM = self.MoriQALY
                self.judge = 0
                self.stage = 6
                self.M_flag = False
            '''





            self.update_transfer()
            #print(self.transfer)



            if week % 52 == 51:
                #3% discount rate

                self.CA1 /= 1.03
                self.CA2 /= 1.03
                self.CB1 /= 1.03
                self.CB2 /= 1.03
                self.CS /= 1.03
                self.CD /= 1.03
                self.CM /= 1.03
                self.QS /= 1.03
                self.UR /= 1.03
                self.UM /= 1.03


                #5% discount rate
                '''
                self.CA /= 1.05
                self.CB /= 1.05
                self.CS /= 1.05
                self.CD /= 1.05
                self.CM /= 1.05
                self.QS /= 1.05
                self.UR /= 1.05
                self.UM /= 1.05
                '''

                '''
                if cumu_T0 + cumu_TA == 0:
                    # failed simu
                    return False
                '''
                cur_annual_rate = (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB) / sum(self.totalstate)
                treatmentfail_annual_rate = cumu_treatment_fail / sum(self.totalstate)
                DST_annual_rate = cumu_DST_in / sum(self.totalstate)
                Annual_case = cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB
                Aresis_annual_rate = (cumu_TA + cumu_PTA) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB + cumu_TC)
                Bresis_annual_rate = (cumu_TB + cumu_PTB) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB + cumu_TC)
                ABresis_annual_rate = (cumu_TAB + cumu_PTAB) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB + cumu_TC)
                Cresis_annual_rate = (cumu_TC) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB + cumu_TC)
                symptomatic_annual_rate = (cumu_T0 + cumu_TA + cumu_TB + cumu_TAB) / (cumu_PT0 + cumu_T0 + cumu_PTA + cumu_TA + cumu_PTB + cumu_TB + cumu_PTAB + cumu_TAB)
                self.annual_case_rate.append(cur_annual_rate)
                self.annual_failtreatment_rate.append(treatmentfail_annual_rate)
                self.annual_DST_rate.append(DST_annual_rate)
                self.annual_incidence_rate.append(cumu_inci / sum(self.totalstate))
                self.annual_resistance_rateA.append(Aresis_annual_rate)
                self.annual_resistance_rateB.append(Bresis_annual_rate)
                self.annual_resistance_rateAB.append(ABresis_annual_rate)
                self.annual_resistance_rateC.append(Cresis_annual_rate)
                self.annual_case.append(Annual_case)

                '''
                if week >= 519 and Aresis_annual_rate ==0:
                    print('error', week)
                '''

                '''
                if Aresis_annual_rate == 0:
                    print('afterAresis_annual_rate: ', week)
                '''

                '''
                if self.stage != 2:
                    self.annual_resistance_rateA.append(Aresis_annual_rate)
                    self.annual_resistance_rateB.append(Bresis_annual_rate)
                    self.annual_resistance_rateAB.append(ABresis_annual_rate)
                    self.annual_resistance_rateC.append(0)
                    self.annual_resistance_rateBC.append(0)
                else:
                    self.annual_resistance_rateA.append(0)
                    self.annual_resistance_rateAB.append(0)
                    self.annual_resistance_rateB.append(Aresis_annual_rate)
                    self.annual_resistance_rateC.append(Cresis_annual_rate)
                    self.annual_resistance_rateBC.append(BCresis_annual_rate)
                '''

                Cost_annual = cumu_annualCost
                QALY_annual = cumu_annualQALY
                self.annul_Cost.append(Cost_annual)
                self.annul_QALY.append(QALY_annual)
                #print('Bresis_annual_rate: ', Bresis_annual_rate)
                #print('ABresis_annual_rate: ', ABresis_annual_rate)

                self.annual_symptomatic_rate.append(symptomatic_annual_rate)

                cumu_PT0 = 0
                cumu_T0 = 0
                cumu_PTA = 0
                cumu_TA = 0
                cumu_PTB = 0
                cumu_TB = 0
                cumu_PTAB = 0
                cumu_TAB = 0
                cumu_inci = 0
                cumu_Asy = 0
                cumu_Sy = 0
                cumu_annualCost = 0
                cumu_annualQALY = 0

                cumu_treatment_fail  = 0
                cumu_DST_in  = 0

        TotalcasesRate = sum(i for i in self.annual_case_rate)
        AveragecaseRate = TotalcasesRate/len(self.annual_case_rate)
        TotalfailtreatmentRate = sum(i for i in self.annual_failtreatment_rate)
        AveragefailtreatmentRate = TotalfailtreatmentRate/len(self.annual_failtreatment_rate)
        TotalDSTRate = sum(i for i in self.annual_DST_rate)
        AverageDSTRate = TotalDSTRate/len(self.annual_DST_rate)
        self.Average_CaseRate.append(AveragecaseRate)
        self.Average_FailtreatmentRate.append(AveragefailtreatmentRate)
        self.Average_DSTRate.append(AverageDSTRate)
        #print('Average_CaseRate: ', self.Average_CaseRate)


        self.Sum_Cost.append(cumu_Totalcost)
        self.Sum_QALY.append(cumu_QALYlost)




        self.Sum_A_cost.append(cumu_A_cost)
        self.Sum_B_cost.append(cumu_B_cost)
        self.Sum_M_cost.append(cumu_M_cost)
        self.Sum_Diagnose_cost.append(cumu_Diagnose_cost)
        self.Sum_Drugtest_cost.append(cumu_Drugtest_cost)
        self.Sum_Sequalae_cost.append(cumu_Sequalae_cost)

        self.Sum_A_qalysloss.append(cumu_A_qalysloss)
        self.Sum_B_qalysloss.append(cumu_B_qalysloss)
        self.Sum_M_qalysloss.append(cumu_M_qalysloss)
        self.Sum_Sequalae_qalysloss.append(cumu_Sequalae_qalysloss)
        self.Sum_Infection_qalyloss.append(cumu_Infection_qalyloss)




        self.FA_Cost.append(cumu_FAcost)
        self.FA_QALY.append(cumu_FAQALY)
        self.FB_Cost.append(cumu_FBcost)
        #self.FB_QALY.append(cumu_BswitchAB)
        '''
        self.FB_QALY.append(cumu_FBQALY)
        self.AB_Cost.append(cumu_ABcost)
        self.AB_QALY.append(cumu_ABQALY)
        '''
        self.AB_Cost.append(cumu_Succ_A)
        self.AB_QALY.append(cumu_Succ_B)
        self.FB_QALY.append(cumu_Succ_M)
        self.Back_number.append(cumu_Succ_back)
        #print('FB_QALY: ', self.FB_QALY)


        self.cumu_N_0.append(cumu_N_0)
        self.cumu_I_0.append(cumu_I_0)
        self.cumu_N_A.append(cumu_N_A)
        self.cumu_I_A.append(cumu_I_A)
        self.cumu_N_B.append(cumu_N_B)
        self.cumu_I_B.append(cumu_I_B)
        self.cumu_N_AB.append(cumu_N_AB)
        self.cumu_I_AB.append(cumu_I_AB)
        self.cumu_PT_0A.append(cumu_PT_0A)
        self.cumu_T_0A.append(cumu_T_0A)
        self.cumu_PT_AA.append(cumu_PT_AA)
        self.cumu_T_AA.append(cumu_T_AA)
        self.cumu_PT_BA.append(cumu_PT_BA)
        self.cumu_T_BA.append(cumu_T_BA)
        self.cumu_PT_ABA.append(cumu_PT_ABA)
        self.cumu_T_ABA.append(cumu_T_ABA)
        self.cumu_F_A.append(cumu_F_A)
        self.cumu_F_ABA.append(cumu_F_ABA)
        self.cumu_PT_0B.append(cumu_PT_0B)
        self.cumu_T_0B.append(cumu_T_0B)
        self.cumu_PT_AB.append(cumu_PT_AB)
        self.cumu_T_AB.append(cumu_T_AB)
        self.cumu_PT_BB.append(cumu_PT_BB)
        self.cumu_T_BB.append(cumu_T_BB)
        self.cumu_PT_ABB.append(cumu_PT_ABB)
        self.cumu_T_ABB.append(cumu_T_ABB)
        self.cumu_F_BB.append(cumu_F_BB)
        self.cumu_F_ABB.append(cumu_F_ABB)
        self.cumu_finalfail.append(cumu_finalfail)

        return True
