"""Span-weighted summaries of one linked chromosome; no tree resampling."""
import numpy as np

def summarize(left,right,values,length,width):
 left,right,values=map(np.asarray,(left,right,values))
 if len(left)==0 or left[0]!=0 or right[-1]!=length or not np.array_equal(left[1:],right[:-1]) or np.any(right<=left):raise ValueError('intervals must tile chromosome')
 edges=np.r_[np.arange(0,length,width),length];mass=np.maximum(0,np.minimum(right[:,None],edges[None,1:])-np.maximum(left[:,None],edges[None,:-1]))
 return (edges[:-1]+edges[1:])/2,values.T@mass/np.diff(edges)
