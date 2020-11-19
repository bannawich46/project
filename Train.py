import os
import gym
import tensorflow
import numpy
import math
import matplotlib
import matplotlib.pyplot as plt

from Env import Watersimulation
from stable_baselines.common.policies import MlpPolicy
from stable_baselines.common import make_vec_env
from stable_baselines import PPO2
from stable_baselines import results_plotter
##log_dir = "bannawich/"
#env = gym.make('Watersimulation')
#env = make_vec_env('Watersimulation', n_envs=4)
#print("debug1")
env = Watersimulation()
#print("debug2")
#env = make_vec_env(lambda: env, n_envs=1)
#env = DummyVecEnv([lambda: env])

#print("debug3")
model = PPO2(MlpPolicy, env, verbose=1)
model.learn(total_timesteps = 60)
##results_plotter.plot_results([log_dir], 60,results_plotter.X_TIMESTEPS, "Graph")
##plt.show()
model.save("ppo2_Watersimulation")

del model
#print("debug4")
model = PPO2.load("ppo2_Watersimulation")
#print("debug5")
observation = [0, 0, 0]

#print("start")
for episode in range(1):
        obs = env.reset()
        for i in range(60):
                #print("Observation :", observation)
                #action, _ = model.action_probability(obs)
                action, _ = model.predict(observation, deterministic = True)
                #action = env.action_space.sample()
                #print("Action :", action)
                #print(action, observation)
                #action, _ = model.predict(obs, env.action_space.low, env.action_space.high)
                #action = numpy.clip(action, env.action_space.low, env.action_space.high)
                observation, reward, dones, info = env.step(action)

                #plt.plot(i,reward)
                #plt.show()

                env.render()

env.close()