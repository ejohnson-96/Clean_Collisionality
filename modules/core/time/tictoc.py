import time

ti = 0  # initial time
def TicTocGenerator():
    # Generator that returns time differences
    tf = time.time()  # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf - ti  # returns the time difference
        tf = ti


TicToc = TicTocGenerator()  # create an instance of the TicTocGen generator


# This will be the main function through which we define both tic() and toc()
def stop_time(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print("\nElapsed time: %f seconds.\n" % tempTimeInterval)


def lap_time(tempBool=True):
    tempTimeInterval = next(TicToc)
    if tempBool:
        print("\nLap time: %f seconds.\n" % tempTimeInterval)


def start_time():
    # Records a time in TicToc, marks the beginning of a time interval
    stop_time(False)


def wait_time(
        wait=5,
):
    if not isinstance(wait, (float, int)):
        raise TypeError(
            "Error: Wait duration in seconds should be of type "
            f"integer or float, instead got type of {type(wait)}"
        )

    time.sleep(wait)
    return

