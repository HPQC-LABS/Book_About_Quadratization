import numpy as np
from scipy.io import loadmat

# hyperparameters
n = 5              # number of variables
aux = True
n_sessions = 100
percentile = 70
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
reset_state.shape = (coeffs_size,)

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

def rhs(s):
    RHS = allbits.dot(s.reshape(coeffs_size,1))
    RHS = RHS - np.min(RHS)
    if aux:
        RHS = np.minimum(RHS[::2],RHS[1::2])    # when using aux
    
    return RHS


def accuracy(s):
    RHS = rhs(s)
    difference = np.abs( np.subtract(RHS,LHS) )
    distance = np.sum(difference)
    preserved_pct = 100*np.size(np.where(difference == 0))/(2*LHS_size)
    
    return preserved_pct, distance


# create agent
from sklearn.neural_network import MLPClassifier
agent = MLPClassifier(
    hidden_layer_sizes=(20,20),
    activation='tanh',
    warm_start=True,  # keep progress between .fit() calls
    max_iter=1,       # make only 1 iteration on each .fit()

)
env = environment()
# initialize agent to the dimension of state and number of actions
agent.fit([env.reset()]*n_actions, range(n_actions))

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
    old_accuracy,temp = accuracy(s)
    
    for t in range(t_max):
        
        # a vector of action probabilities in current state
        probs = agent.predict_proba([s])[0]
        
        #choose 100 of those
        a = np.random.choice(n_actions, p = probs)
        
        new_s = env.step(a)
        
        new_accuracy,dist = accuracy(new_s)
        
        # Reward System
        r = (new_accuracy - old_accuracy)   # reward/punishment for number of states preserved        
        if dist > LHS_size:
            r -= (dist - LHS_size)          # punishment for distance away from LHS
        if np.max(np.abs(new_s)) >= 4:
            r -= 10                         # punishment for large coefficients
        
        if new_accuracy >= 80:
            flag = True
            for k in range(len(good)):
                if (good[k] == new_s).all:
                    flag = False
            if flag:
                good.append(new_s)
            print(new_s,new_accuracy)
            #print(rhs(new_s))
            #wait = input("press enter to continue")
        
        states.append(s)
        actions.append(a)
        total_reward += r
        
        s = new_s
        old_accuracy = new_accuracy
        
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
    mask = rewards_batch > reward_threshold
    
    elite_states  = [ item for i in range(len(states_batch))  for item in states_batch[i]  if mask[i] ]
    elite_actions = [ item for i in range(len(actions_batch)) for item in actions_batch[i] if mask[i] ]
    
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
for i in range(200):
    #generate new sessions
    sessions = [ generate_session() for _ in range(n_sessions) ]

    batch_states,batch_actions,batch_rewards = map(np.array, zip(*sessions))
    #print(batch_rewards)
    #wait = input("rewards")
    elite_states, elite_actions = select_elites(batch_states, batch_actions, batch_rewards, percentile)

    agent.fit(elite_states, elite_actions)

    show_progress(batch_rewards, percentile, reward_range=[0,np.max(batch_rewards)])

#    if np.mean(batch_rewards)> 50:
 #       print("You Win!!! You may stop training now")

