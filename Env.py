import numpy as np
import zmq
import pickle
import gym
import math
import matplotlib.pyplot as plt
from gym import spaces

context = zmq.Context()
socket = context.socket(zmq.REQ)
socket.connect("tcp://localhost:7256")


class Watersimulation(gym.Env):
    def __init__(self):
        super(Watersimulation, self).__init__()
        self.end = 0
        self.action_space = gym.spaces.Box(low = np.array([0,0]), high = np.array([16,16]), dtype = np.int64)
        self.observation_space = gym.spaces.Box(low = np.array([0.0, 0.0, 0.0]), high = np.array([100.0, 100.0, 100.0]), dtype = np.float64)

    def step(self, action):
        self.end = self.end + 1
        dones = bool(self.end == 1)
        #print(type(action))
        #action1 = self.action_space.sample()[1]
        #action2 = self.action_space.sample()[-1]
        #action = [action1, action2]
        #print(action)
        #actiongo = [action, 0]
        print("action",action, type(action))
        info = {}
        action = [0, 0]
        r_num = random... # 0-1
        if r_num <= 0.1:
            action = self.action_space.sample()

        go1 = pickle.dumps(action, protocol = 2)
        socket.send(go1)
        
        
        come1 = socket.recv()
        informationcome1 = pickle.loads(come1)
        [rewards, check1, check2, check3] = informationcome1
        observation = [check1, check2, check3]
        reward = rewards
        #print(check3,",")
        #print(reward,",")
        return observation, reward, dones, info
    def reset(self):
        pass

    def render(self, mode='human'):
        pass

    def close(self):
        informationgo2 = "close"
        go2 = pickle.dumps(informationgo2, protocol = 2)
        socket.send(go2)