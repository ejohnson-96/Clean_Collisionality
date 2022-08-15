import time
import datetime
from modules.core.variables.num_man import valid_int

def validate_time(
        alarm_time,
):
    if len(alarm_time) != 11:
        return "Invalid time format! Format HH:MM:SS is required," \
               f" instead got {alarm_time}. Please try again..."
    else:
        if int(alarm_time[0:2]) > 24:
            return "Invalid HOUR format! Hour argument provided " \
                   f"{alarm_time[0:2]}. Please try again..."
        elif int(alarm_time[3:5]) > 59:
            return "Invalid MINUTE format! Minute argument provided " \
                   f"{alarm_time[3:5]}. Please try again..."
        elif int(alarm_time[6:8]) > 59:
            return "Invalid SECOND format! Second argument provided " \
                   f"{alarm_time[6:8]}. Please try again..."
        else:
            return alarm_time

def countdown(
        t_,
):
    t_ = valid_int(t_)
    if t_ < 0:
        raise ValueError(
            "Error: Countdown number cannot be negative."
        )
    elif t_ == 0:
        raise ValueError(
            "Warning: Countdown number set to zero."
        )
    else:
        while t_:
            mins, secs = divmod(t_, 60)
            timer = '{:02d}:{:02d}'.format(mins, secs)
            print(timer, end="\r")
            time.sleep(1)
            t_ -= 1
    return


def alarm_clock(
        alarm_time,
):
    validate_time(alarm_time.lower())
    print(f"Setting alarm for {alarm_time}...")

    alarm_hour = alarm_time[0:2]
    alarm_min = alarm_time[3:5]
    alarm_sec = alarm_time[6:8]
    alarm_period = alarm_time[9:].upper()

    now = datetime.now()

    h = 0
    while h < 1:
        if alarm_period == now.strftime("%p"):
            if alarm_hour == now.strftime("%I"):
                if alarm_min == now.strftime("%M"):
                    if alarm_sec == now.strftime("%S"):
                        print("Wake Up!")
                        h = 1
    return