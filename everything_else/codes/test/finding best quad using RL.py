import numpy as np
from scipy.io import loadmat

# hyperparameters
n = 7              # number of variables
aux = True
conflict_threshold = 0.1
n_sessions = 100
t_max = 50
c_conflict = 1
c_distance = 1
c_entropy = 0.4
percentile = 70
hidden_layers = (20,30,40)
np.random.seed(1)
# -----------------------------------------------------------------------------------------

coeffs_size = int( n*(n+1)/2 )
n_actions = 2*coeffs_size + 1     # number of actions

if aux:
    LHS_size = 2**(n-1)
else:
    LHS_size = 2**n

LHS = loadmat('data.mat')['LHS']
LHS.shape = (LHS_size,1)

allbits = loadmat('data.mat')['allbits']
allbits = np.transpose(allbits)

reset_state = loadmat('data.mat')['reset_state']
reset_state.shape = (coeffs_size,-1)
reset_state = np.repeat(reset_state,n_sessions,1)

one_hot = np.zeros((coeffs_size,n_actions))
for i in range(coeffs_size):
    one_hot[i,i] = 1
    one_hot[i,coeffs_size+i] = -1

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
    #if np.max(np.abs(new_s)) >= 4:
    #    r -= 1                           # punishment for large coefficients
    
    return c_conflict * r_conflict + c_distance * r_distance + c_entropy * r_entropy

def rhs(s):
    RHS = allbits.dot(s)
    RHS = RHS - np.min(RHS,0)
    if aux:
        RHS = np.minimum(RHS[::2],RHS[1::2])    # when using aux
    
    return RHS


def accuracy(s):
    RHS = rhs(s)
    difference = np.abs( np.subtract(RHS,LHS) )
    distance = np.sum(difference,0)
    distance = np.minimum(1,distance/12) #LHS_size? or maybe changing at runtime?
    conflict = np.sum(difference != 0,0) / LHS_size
    
    return conflict, distance


# create agent
from sklearn.neural_network import MLPClassifier
agent = MLPClassifier(
    hidden_layer_sizes = hidden_layers,
    activation='tanh',
    warm_start=False,#True,  # keep progress between .fit() calls
    max_iter=1,       # make only 1 iteration on each .fit()

)
env = environment()
# initialize agent to the dimension of state and number of actions
init_training = loadmat('data.mat')['input']
size = init_training.shape[0]

train_data = np.subtract(init_training,one_hot[:,0])
target = np.ones((size,1))*0
for action in range(1,n_actions):
    train_data = np.append( train_data, np.subtract(init_training,one_hot[:,action]) )
    target = np.append( target, np.ones((size,1))*action )

train_data = train_data.reshape(-1,coeffs_size)

agent.fit(train_data, target)


def generate_session(t_max = 50):
    states,actions = [],[]
    total_reward = 0
   
    s = env.reset()
        
    for t in range(t_max):
        
        probs = agent.predict_proba(np.transpose(s))
        
        # choose actions w.r.t. probs
        c = probs.cumsum(axis = 1)
        u = np.random.rand(len(c),1)
        a = (u < c).argmax(axis = 1)
        
        new_s = env.step(a)
        
        conflict,dist = accuracy(new_s)
        entropy = -np.sum( np.multiply(probs,np.log(probs)/np.log(coeffs_size)), 1)
        
        r = reward(conflict, dist, entropy)
        
        if any(conflict <= conflict_threshold):
            #flag = True
            #for k in range(len(good)):
            #    if (good[k] == new_s).all:
            #        flag = False
            #if flag:
            #   good.append(new_s)
            #print(conflict)
            print( new_s[:,conflict <= conflict_threshold], 100*(1-conflict[conflict <= conflict_threshold]) )
            wait = input("press enter to continue")
        
        states.append(np.transpose(s))
        actions.append(a)
        total_reward += r
        
        s = new_s
        
    return states, actions, total_reward


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


from IPython.display import clear_output

def show_progress(batch_rewards, percentile, reward_range=[0,100]):
    
    mean_reward, threshold = np.mean(batch_rewards), np.percentile(batch_rewards, percentile)
    Log.append([int(mean_reward), int(threshold)])

    clear_output(True)
    print("mean reward = %.3f, threshold=%.3f"%(mean_reward, threshold))
    """
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
    """

Log  = []
good = []
for i in range(500):
    #sessions = [ generate_session() for _ in range(n_sessions) ]
    sessions = generate_session(t_max)
    
    #batch_states,batch_actions,batch_rewards = map(np.array, zip(*sessions))
    batch_states,batch_actions,batch_rewards = sessions[0], sessions[1], sessions[2]
    #print(batch_rewards)
    #wait = input("rewards")
    elite_states, elite_actions = select_elites(batch_states, batch_actions, batch_rewards, percentile)
    
    agent.partial_fit(elite_states, elite_actions)
    
    show_progress(batch_rewards, percentile, reward_range=[0,np.max(batch_rewards)])
    
#    if np.mean(batch_rewards)> 50:
 #       print("You Win!!! You may stop training now")

