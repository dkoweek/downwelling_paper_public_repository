import pandas as pd
import numpy as np
import numpy.matlib
from scipy import stats
from sklearn import linear_model
import statsmodels.api as sm
from statsmodels.tsa.api import ExponentialSmoothing, SimpleExpSmoothing, Holt








def bin_sets(x,y,x_step):

  x_b = np.arange(np.nanmin(x),np.nanmax(x)+np.nanmin(x)/1.e6,x_step)
  y_o = np.zeros(len(x_b)-1)
  x_o = np.zeros(len(x_b)-1)
  y_s = np.zeros(len(x_b)-1)
  for i in range(0,len(x_b)-1):
    subset = np.logical_and(x >= x_b[i],x<x_b[i+1])
    y_o[i] = np.nanmean(y[subset])
    y_s[i] = np.nanstd(y[subset])
    x_o[i] = np.nanmean(x[subset])

  return x_o,y_o,y_s



df = pd.read_csv('data/oxygenation_expt_sensor_data.csv',low_memory=False, converters={'serial_number': lambda x: str(x)})
print('loaded')

sensors = np.array(df['serial_number'])
mooring = np.array(df['Mooring']).astype(str)
un_sensors = np.unique(sensors)
print(un_sensors)
sensor_radius = np.array([np.array(df['Radius'][val == sensors])[0] for val in un_sensors])
sensor_mab = np.array([np.array(df['Meters_above_bottom'][val == sensors])[0] for val in un_sensors])
radius = np.array(df['Radius'])
angle = np.array(df['Angle_CW_from_N'])


oxygen = np.array(df['O2'])
meters_above_bottom = np.array(df['Meters_above_bottom'])
print('grabbed components')

time = pd.to_datetime(df['datetime']) - pd.Timedelta(hours=7)
print('datetime components')
minute = time.dt.hour*60+time.dt.minute
day = np.array(df['treatment_ID_local'],dtype='str')

un_day = np.unique(day)
print(un_day)

pumping = (np.array(df['pumps_on']) == 1)
would_have_pumping = (np.array(df['pumps_on']) == 0)


un_sensor_subset = np.unique(sensors[np.logical_or.reduce((radius == 2, radius == 5, radius == 8))])

agg_window = 10 


def fill_in_day(dat,length):
  tmp_ = np.zeros(length)
  tmp_[:] = np.nan
  tmp_[-len(dat):] = dat
  return tmp_
 


np.random.seed(13)

for _time in [30, 60, 120]:

    pump_buffer = int(_time / agg_window)

    all_deviations = np.zeros((len(un_sensor_subset),len(un_day),4))
    all_pumping = np.zeros((len(un_sensor_subset),len(un_day),4))
    
   
    for _n in range(len(un_sensor_subset)):
    
      day_ox = []
      day_pump = []
      day_would_have_pump = []
      day_min = []
      store_day = []
    
      is_pump_on = []
    
      subset = np.logical_and.reduce((sensors == un_sensor_subset[_n],\
                                      minute > 8 * 60))
    
      l_un_day = np.unique(day[subset])
    
      length = []
      for _d in range(0,len(l_un_day)):
        day_subset = np.logical_and(subset,day == l_un_day[_d])
    
        mean_min,mean_ox,std_ox = bin_sets(minute[day_subset],oxygen[day_subset],agg_window)
        mean_min,mean_pump,std_pump = bin_sets(minute[day_subset],pumping[day_subset].astype(float),agg_window)
        mean_min,mean_would_have_pump,std_would_have_pump = bin_sets(minute[day_subset],would_have_pumping[day_subset].astype(float),agg_window)
    
        day_ox.append(mean_ox)
        day_pump.append(mean_pump)
        day_would_have_pump.append(mean_would_have_pump)
        day_min.append(mean_min)
        length.append(len(day_pump[-1]))
        store_day.append(l_un_day[_d])
    
      mode_len = stats.mode(length)[0][0]
      for _d in range(len(length)-1,-1,-1):
        if (length[_d] != mode_len and l_un_day[_d] != 'C1_PDT'):
          day_ox.pop(_d)
          day_pump.pop(_d)
          day_would_have_pump.pop(_d)
          length.pop(_d)
          day_min.pop(_d)
          store_day.pop(_d)
        elif (l_un_day[_d] == 'C1_PDT'): #we know that day 0 didn't have data untill 1700 local
          day_ox[_d] = fill_in_day(day_ox[_d],mode_len)
          day_pump[_d] = fill_in_day(day_pump[_d],mode_len)
          day_would_have_pump[_d] = fill_in_day(day_would_have_pump[_d],mode_len)
          day_min[_d] = fill_in_day(day_min[_d],mode_len)
          print('not popping C1_PDT')
    
      day_ox = np.vstack(day_ox)
      day_pump = np.vstack(day_pump)
      day_would_have_pump = np.vstack(day_would_have_pump)
      day_min = np.vstack(day_min)
      store_day = np.array(store_day)
    
    
      pred_ox = np.zeros(day_ox.shape)
      pred_ox[:] = np.nan
      day_has_pumping = np.zeros(pred_ox.shape[0]).astype(bool)
    
      for _d in range(day_ox.shape[0]):
    
        if (np.nansum(day_pump[_d]) == 0):
          first_pump = np.nanargmax(day_would_have_pump[_d])
          last_pump = day_min.shape[1]-np.nanargmax(day_would_have_pump[_d][::-1])
          day_has_pumping[_d] = False
        else:
          first_pump = np.nanargmax(day_pump[_d])
          last_pump = day_min.shape[1]-np.nanargmax(day_pump[_d][::-1])
          day_has_pumping[_d] = True
        pre_pump_mean = np.nanmean(day_ox[_d,first_pump-pump_buffer:first_pump])
        post_pump_mean = np.nanmean(day_ox[_d,last_pump+1:last_pump+1+pump_buffer])
        pred_len = last_pump-first_pump
    
        pred_ox[_d,first_pump:first_pump+pred_len] = np.nanmean(day_ox[_d,first_pump-pump_buffer:first_pump-1]) + (post_pump_mean - pre_pump_mean)/float(last_pump-first_pump+pump_buffer) * (np.arange(pred_len) + pump_buffer/2.)
    
      deviation = day_ox - pred_ox
    
      for _d in range(len(store_day)):
    
        first_pump = np.argmax(np.isnan(deviation[_d,:]) == False)
        last_pump = deviation.shape[1] - np.argmax((np.isnan(deviation[_d,:]) == False)[::-1])+1
        time_diff = last_pump-first_pump
    
        _i1 = un_sensor_subset == un_sensor_subset[_n]
        _i2 = un_day == store_day[_d]
        all_deviations[_i1,_i2,0] = np.nanmean(deviation[_d,first_pump:last_pump])
        all_pumping[_i1,_i2,0] = np.nanmean(day_pump[_d,first_pump:last_pump])
    
        all_deviations[_i1,_i2,1] = np.nanmean(deviation[_d,first_pump:first_pump + int(time_diff/3.)])
        all_pumping[_i1,_i2,1] = np.nanmean(day_pump[_d,first_pump:first_pump + int(time_diff/3.)])
    
        all_deviations[_i1,_i2,2] = np.nanmean(deviation[_d,first_pump + int(time_diff/3.):first_pump + int(time_diff*2/3.)])
        all_pumping[_i1,_i2,2] = np.nanmean(day_pump[_d,first_pump + int(time_diff/3.):first_pump + int(time_diff*2/3.)])
    
        all_deviations[_i1,_i2,3] = np.nanmean(deviation[_d,first_pump + int(time_diff*2/3.):last_pump])
        all_pumping[_i1,_i2,3] = np.nanmean(day_pump[_d,first_pump + int(time_diff*2/3.):last_pump])
    
    #np.savez('munged/deviation_dat.npz',all_deviations=all_deviations,all_pumping=all_pumping,day=l_un_day,sensor=un_sensor_subset)    

    # Convert to output
    all_deviation[all_deviation == 0] = np.nan

    sensor_angle = []
    sensor_radius = []
    x_pos = []
    y_pos = []
    mab = []
    
    for _n in range(len(un_sensor_subset)):
      ind = np.argmax(sensors == un_sensor_subset[_n])
      sensor_angle.append(angle[ind])
      sensor_radius.append(radius[ind])
      mab.append(meters_above_bottom[ind])
    
      x_pos.append(np.sin(sensor_angle[-1]*np.pi/180.)*sensor_radius[-1])
      y_pos.append(np.cos(sensor_angle[-1]*np.pi/180.)*sensor_radius[-1])
    
      print((un_sensor_subset[_n],sensor_radius[-1],sensor_angle[-1],x_pos[-1],y_pos[-1],np.nanmean(all_deviation[_n,:])))
    
    x_pos = np.array(x_pos)
    y_pos = np.array(y_pos)
    mab = np.array(mab)
    print((len(x_pos),len(y_pos),len(mab)))
    
    
    df = pd.DataFrame(data=un_sensor_subset,columns=['Sensor'])
    df['x_pos'] = x_pos
    df['y_pos'] = y_pos
    df['mab'] = mab
    dev_day = l_un_day
    for i in range(all_deviation.shape[1]):
      df['Day ' + str(dev_day[i]) + ' dev.'] = all_deviation[:,i,0]
    for i in range(all_deviation.shape[1]):
      df['Day ' + str(dev_day[i]) + ' dev. hour 1'] = all_deviation[:,i,1]
    for i in range(all_deviation.shape[1]):
      df['Day ' + str(dev_day[i]) + ' dev. hour 2'] = all_deviation[:,i,2]
    for i in range(all_deviation.shape[1]):
      df['Day ' + str(dev_day[i]) + ' dev. hour 3'] = all_deviation[:,i,3]
    for i in range(all_pumping.shape[1]):
      df['Day ' + str(dev_day[i]) + ' pumping'] = all_pumping[:,i,0]
    for i in range(all_pumping.shape[1]):
      df['Day ' + str(dev_day[i]) + ' pumping hour 1'] = all_pumping[:,i,1]
    for i in range(all_pumping.shape[1]):
      df['Day ' + str(dev_day[i]) + ' pumping hour 2'] = all_pumping[:,i,2]
    for i in range(all_pumping.shape[1]):
      df['Day ' + str(dev_day[i]) + ' pumping hour 3'] = all_pumping[:,i,3]
    
    if (_time == 30):
        time_tag = 'half_hour'
    elif (_time == 60):
        time_tag = '1h'
    elif (_time == 120):
        time_tag = '2h'
    else:
        time_tag = str(_time) + 'm'
    df.to_csv('data/mapped_pumping_data_' + time_tag + '.csv',index=False,sep=',')
     

















