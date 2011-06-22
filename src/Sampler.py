from abc import ABCMeta, abstractmethod

class Sampler:
	__metaclass__ = ABCMeta

	@abstractmethod
	def subsample(self, se):
		return NotImplemented

