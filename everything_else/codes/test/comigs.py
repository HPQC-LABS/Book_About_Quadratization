import numpy as np
from scipy.io import savemat # ,loadmat
import matlab.engine
from sklearn.neural_network import MLPClassifier

# hyperparameters
#n = 7              # number of variables including aux
#aux = 1
fileID = "LHS.txt"
num_cases = 1
new_input = True

n_sessions = 100
t_max = 150

c_conflict = 1
c_distance = 1
c_entropy  = 0.25

c_num_of_terms = 1
c_strength = 2
c_accuracy = 10000

percentile = 70
max_comigs = 15
comigs_each_run = 1
# -----------------------------------------------------------------------------------------

# functions/classes

class environment:
    
    def reset(self):
        self.state = reset_state
        return self.state
    
    def step(self, action):
        self.state = np.add( self.state, one_hot[:,action])
        return self.state

def reward(conflict, dist, entropy):
    r_conflict = 1 - conflict**(1/3)             # "preserving states" reward
    r_distance = 1 - dist**(1/3)                 # "staying close"     reward
    r_entropy  = entropy                         # "being an explorer" reward
    
    return c_conflict * r_conflict + c_distance * r_distance + c_entropy * r_entropy

def rhs(s):
    RHS = allbits.dot(s)
    RHS = RHS - np.min(RHS,0) + min(LHS)
    if aux == 1:
        RHS  = np.minimum(RHS[::2],RHS[1::2])
    if aux == 2:
        min1 = np.minimum(RHS[::4],RHS[1::4])
        min2 = np.minimum(RHS[2::4],RHS[3::4])
        RHS  = np.minimum(min1,min2)
    
    return RHS

def accuracy(s):
    
    RHS = rhs(s)
    difference = np.subtract(RHS,LHS)
    flag = (np.sum(difference < 0,0) != 0)
    
    conflict = np.sum(difference[mask] != 0,0) / sum(mask)
    
    difference = difference[difference < 0]
    difference = np.abs(difference)
    distance = np.multiply( np.sum(difference,0) ,flag.reshape(-1,) )
    distance = np.minimum(1, distance/sum(mask)) #LHS_size? or maybe changing at runtime?
    
    return conflict, distance, flag

def generate_session(t_max = 50):
    global count, worst_score, worst_index, stop_condition, good, found
    states,actions = [],[]
    total_reward = 0
    mini = 100
    #r_best = -1000 * np.ones((n_sessions,))
    s = env.reset()
    #conflict,dist = accuracy(s)
    #old_reward = reward(conflict, dist, 0)
    
    for t in range(t_max):
        
        probs = agent.predict_proba(np.transpose(s))
        
        # choose actions w.r.t. probs
        c = probs.cumsum(axis = 1)
        u = np.random.rand(len(c),1)
        a = (u < c).argmax(axis = 1)
        
        new_s = env.step(a)
        
        conflict, dist, flag = accuracy(new_s)
        
        entropy = -np.sum( np.multiply(probs,np.log(probs)/np.log(coeffs_size)), 1)
        
        #new_reward = reward(conflict, dist, entropy)
        #r = new_reward - old_reward
        #old_reward = new_reward
        
        r = reward(conflict, dist, entropy)
        r[flag] -= 10
        
        condition = ( flag == 0 )
        
        if any(condition):
            coeffs_strength = np.sum(np.abs(new_s[:,condition]),0)
            coeffs_number = np.sum( new_s[:,condition] != 0 , 0 )
            score = c_accuracy*(1 - conflict[condition]) - c_strength*coeffs_strength - c_num_of_terms*coeffs_number
            score_condition = ( score > worst_score )
            
            if any(score_condition):
                stop_condition = 0
                found = 1
                #flag = True
                #for k in range(len(good)):
                #    if (good[k] == new_s[:,temp]).all():
                #        flag = False
                #if flag:
                #   good.append(new_s[:,temp])
                #print(rhs(new_s[:,temp]))
                coeffs = np.transpose(new_s[:,condition])
                coeffs = coeffs[score_condition]
                score = score[score_condition]
                #print(j,max(score),worst_score)
                for i in range(coeffs.shape[0]):
                    if score[i] > worst_score:
                        if count < comigs_each_run:
                            good[num_COMIGs,count] = coeffs[i]
                            count = count + 1
                        else:
                            good[num_COMIGs,worst_index] = coeffs[i]
                        
                        if count == comigs_each_run:
                            worst_score = 10000
                            worst_index = 10000
                            for it in range(comigs_each_run):
                                coef = good[num_COMIGs,it]
                                conflict_current, _, _ = accuracy(coef.reshape(-1,1))
                                coeffs_strength = np.sum(np.abs(coef))
                                coeffs_number = np.sum(coef != 0)
                                score_current = c_accuracy*(1 - conflict_current) - c_strength*coeffs_strength - c_num_of_terms*coeffs_number
                                if score_current < worst_score:
                                    worst_score = score_current
                                    worst_index = it
        
        if min(conflict) < mini:
            mini = min(conflict)
        
        states.append(np.transpose(s))
        actions.append(a)
        total_reward += r
        
        s = new_s
    
    stop_condition += 1
    flag = False
    if stop_condition == 40:
        flag = True
    #print(1-mini)
    return states, actions, total_reward, flag, found

def select_elites(states_batch,actions_batch,rewards_batch,percentile=50):
    """
    Select states and actions from games that have rewards >= percentile
    :param states_batch: list of lists of states, states_batch[session_i][t]
    :param actions_batch: list of lists of actions, actions_batch[session_i][t]
    :param rewards_batch: list of rewards, rewards_batch[session_i][t]
    
    :returns: elite_states,elite_actions, both 1D lists of states and respective actions from elite sessions
    """
    
    reward_threshold = np.percentile(rewards_batch, percentile)
    mask = (rewards_batch > reward_threshold)
    
    states_batch = np.array(states_batch).transpose(1,0,2).reshape(n_sessions,-1)
    elite_states = states_batch[mask].reshape(-1,coeffs_size)
    
    actions_batch = np.transpose(actions_batch)
    elite_actions = actions_batch[mask].reshape(-1,1)
    
    return elite_states, elite_actions

def show_progress(batch_rewards, percentile):
    
    mean_reward, threshold = np.mean(batch_rewards), np.percentile(batch_rewards, percentile)
    #print("mean reward = %.3f, threshold=%.3f"%(mean_reward, threshold))
    #Log.append([int(mean_reward), int(threshold)])
    
    print("mean reward = %.3f, threshold=%.3f"%(mean_reward, threshold))
    return mean_reward
    '''
    plt.figure(figsize=[8,4])
    plt.subplot(1,2,1)
    plt.plot(list(zip(*log))[0], label='Mean rewards')
    plt.plot(list(zip(*log))[1], label='Reward thresholds')
    plt.legend()
    plt.grid()
    
    plt.subplot(1,2,2)
    plt.hist(batch_rewards, range=reward_range);
    plt.vlines([np.percentile(batch_rewards, percentile)], [0], [100], label="percentile", color='red')
    plt.legend()
    plt.grid()

    plt.show()
    '''
#-----------------------------------------------------------------------------------------------

# Run for multiple LHSs

for Case in range(num_cases):
    # Read new input
    
    if new_input:
        
        file1 = open(fileID,"r")
        
        LHS_string = file1.readline()
        aux = int( file1.readline() )
        file1.close()
        
        # Run Matlab to pretrain RL model
        eng = matlab.engine.start_matlab()
        #eng.brute_force_for_RL_input_data(nargout = 0)
        [init_training, LHS, allbits, reset_state, n] = eng.pretrain(aux, LHS_string, nargout = 5)
        
        coeffs_size = int( n*(n+1)/2 )
        n_actions = 2*coeffs_size + 1   # number of actions
        
        LHS_size = 2**(n-aux)
        
        #LHS = loadmat('data.mat')['LHS']
        LHS = np.array(LHS)
        LHS.shape = (LHS_size,1)
        
        #allbits = loadmat('data.mat')['allbits']
        allbits = np.array(allbits)
        allbits = np.transpose(allbits)
        
        #reset_state = loadmat('data.mat')['reset_state']
        reset_state = np.array(reset_state)
        reset_state.shape = (coeffs_size,-1)
        reset_state = np.repeat(reset_state,n_sessions,1)
        
        one_hot = np.zeros((coeffs_size,n_actions))
        for i in range(coeffs_size):
            one_hot[i,i] = 1
            one_hot[i,coeffs_size+i] = -1
    
    # create agent
    
    agent = MLPClassifier(
        hidden_layer_sizes = (coeffs_size,2*coeffs_size),
        activation='tanh',
        warm_start=False,#True,  # keep progress between .fit() calls
        max_iter=1,       # make only 1 iteration on each .fit()
    )
    
    env = environment()
    
    #init_training = loadmat('data.mat')['input']
    init_training = np.array(init_training)
    size = init_training.shape[0]
    
    train_data = np.subtract(init_training,one_hot[:,0])
    target = np.ones((size,1))*0
    for action in range(1,n_actions):
        train_data = np.append( train_data, np.subtract(init_training,one_hot[:,action]) )
        target = np.append( target, np.ones((size,1))*action )
    
    train_data = train_data.reshape(-1,coeffs_size)
    
    #------------------------------------------------------------------------------------------------
    # Finding COMIGs
    
    best_num_of_comigs = 1000
    for Run in range(3):
        good = np.zeros((max_comigs,comigs_each_run,coeffs_size))
        good_ind = []
        score = 0
        
        for num_COMIGs in range(max_comigs):
            
            coef = np.zeros((max_comigs,coeffs_size))
            for i in range(num_COMIGs):
                coef[i] = good[i,good_ind[i]]
            
            mask = np.array([True]*LHS_size)
            mask = mask.reshape(-1,1)
            mask_preserve = np.zeros((max_comigs,LHS_size,1))
            
            for i in np.arange(num_COMIGs):
                RHS = rhs(coef[i].reshape(-1,)).reshape(-1,1)
                mask_preserve[i] = np.subtract(RHS,LHS) != 0
                mask = np.multiply(mask, mask_preserve[i])
            
            mask = mask.reshape(-1,)
            mask = mask == 1
            if sum(mask) == 0:
                break
            print('states left = ', sum(mask))
            
            count = 0
            worst_score = -10000
            worst_index = -10000
            stop_counter = 0
            counter = 0
            found = 0
            while counter < 4:
                stop_condition = 0
                np.random.seed()
                agent.fit(train_data,target) # initialize agent
                for j in range(100):
                    #mean_reward_max_in = -10000
                    
                    sessions = generate_session(t_max)
                    
                    batch_states,batch_actions,batch_rewards, flag = sessions[0], sessions[1], sessions[2], sessions[3]
                    
                    if flag:
                        break
                    
                    elite_states, elite_actions = select_elites(batch_states, batch_actions, batch_rewards, percentile)
                    
                    agent.partial_fit(elite_states, elite_actions)
                    
                    '''
                    flag = False
                    if j % 10 == 0:
                        flag = True
                    
                    mean_reward = show_progress(batch_rewards, percentile)#, flag)
                    
                    if mean_reward > mean_reward_max_in:
                        mean_reward_max_in = mean_reward
                
                if (mean_reward_max_in / (c_distance + c_entropy) / t_max) > mean_reward_max:
                    mean_reward_max = mean_reward_max_in
                    i_max = [t_max, c_distance, c_entropy]
                    print(i_max)
                        
                    if np.mean(batch_rewards)> 50:
                        print("You Win!!! You may stop training now")
                '''
                if counter == 3 and found == 0 and stop_counter < 6:
                    counter = 2
                    stop_counter += 1
                counter = counter + 1
            
            print('Case: ',Case,'Run: ',Run,'Comig: ',num_COMIGs)
            best_score = -10000
            best_index = -10000
            for it in range(comigs_each_run):
                coef_temp = good[num_COMIGs,it]
                conflict_current, _, _ = accuracy(coef_temp.reshape(-1,1))
                coeffs_strength = np.sum(np.abs(coef_temp))
                coeffs_number = np.sum(coef_temp != 0)
                score_current = c_accuracy*(1 - conflict_current) - c_strength*coeffs_strength - c_num_of_terms*coeffs_number
                if score_current > best_score:
                    best_score = score_current
                    best_index = it
            good_ind.append(best_index)
            score = score + best_score
        print('num of comigs = ',num_COMIGs)
        if num_COMIGs < best_num_of_comigs:
            best_num_of_comigs = num_COMIGs
            Best_score = score
            best_coef = coef[:num_COMIGs]
            savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
        elif num_COMIGs == best_num_of_comigs:
            if score > Best_score:                
                best_coef = coef[:num_COMIGs]
                savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
    
    print('Done!')
    print('best num of comigs = ',best_num_of_comigs,'\n')
    
    # verify coef and get const_terms
    [verified, const_terms] = eng.verify(n, aux, LHS.tolist(), allbits.tolist(), best_coef.tolist(), nargout = 2)
    
    file1 = open(fileID,"a")
    
    if (np.min(rhs(np.transpose(best_coef)),1) != np.transpose(LHS)).any():
        file1.write('\nNot verified!\n')
    else:
        file1.write('\nVerified!\n')
    
    for m in range(best_coef.shape[0]):
        vec = best_coef[m,:]
        k = 0
        for i in range(1,n):
            for j in range(i+1,n+1):
                if vec[k] != 0:
                    if vec[k] == 1:
                        file1.write(' + b_{%d}b_{%d}' % (i, j) )
                    elif vec[k] == -1:
                        file1.write(' - b_{%d}b_{%d}' % (i, j) )
                    else:
                        file1.write(' %+d b_{%d}b_{%d}' % (vec[k], i, j) )
                k += 1
    
        for i in range(1,n+1):
            if vec[k] != 0:
                if vec[k] == 1:
                    file1.write(' + b_{%d}' % i)
                elif vec[k] == -1:
                    file1.write(' - b_{%d}' % i)
                else:
                    file1.write(' %+d b_{%d}' % (vec[k], i) )
            k += 1
        
        if np.floor(const_terms[m]) != 0:
            file1.write(' %+d' % np.floor(const_terms[m]) )
        file1.write('\n')
    
    file1.close()

