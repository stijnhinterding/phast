/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#include "fakecommunicator.h"

#include <gsl/gsl_randist.h>
#include <math.h>

#include <iostream>

#if defined(_WIN32) && !defined(NDEBUG)
// For GetAsyncKeystate
#include <windows.h>
#endif

FakeCommunicator::~FakeCommunicator()
{
    gsl_rng_free(this->r);
}

FakeCommunicator::FakeCommunicator(double timeunit_seconds, double decayrate_inv_seconds, double pulse_period_seconds, double detection_probability, uint64_t sync_divider,
                                   int64_t chan1_induced_delay, int64_t chan2_induced_delay, int64_t n_emitters, double sim_speedup_factor) :
    r(nullptr),
    timeunit_seconds(timeunit_seconds),
    decayrate_inv_seconds(decayrate_inv_seconds),
    pulse_period_seconds(pulse_period_seconds),
    detection_probability(detection_probability),
    prev_time(std::chrono::high_resolution_clock::now()),
    first_run(true),
    last_largest_time(0),
    last_pulse_time(0),
    pulse_period_timeunits(this->pulse_period_seconds/this->timeunit_seconds),
    last_pulse_published(0),
    sync_counter(0),
    sync_divider(sync_divider),
    chan1_induced_delay(chan1_induced_delay),
    chan2_induced_delay(chan2_induced_delay),
    descriptor("fake-device"),
    chan0_enabled(true),
    chan1_enabled(true),
    chan2_enabled(true),
    n_emitters(n_emitters),
    sim_speedup_factor(sim_speedup_factor),
#if defined(_WIN32) && !defined(NDEBUG)
    is_running(false)
#else
    is_running(true)
#endif
{
    this->r = gsl_rng_alloc (gsl_rng_ranlxs2);
    gsl_rng_set(r, 0);

    this->n_events[0] = 0;
    this->n_events[1] = 0;
    this->n_events[2] = 0;

    this->first_times[0] = 0;
    this->first_times[1] = 0;
    this->first_times[2] = 0;

    this->last_times[0] = 0;
    this->last_times[1] = 0;
    this->last_times[2] = 0;
}

void FakeCommunicator::DisableChan(uint64_t ID)
{
    if (ID == 0) {
        this->chan0_enabled = false;
    } else if (ID == 1) {
        this->chan1_enabled = false;
    } else if (ID == 2) {
        this->chan2_enabled = false;
    }
}

void FakeCommunicator::gen_pulses(std::vector<int64_t>* timestamps, std::vector<uint8_t>* channel_IDs, int64_t duration)
{
    std::lock_guard<std::mutex> l(this->m);

#if defined(_WIN32) && !defined(NDEBUG)
    if (!this->is_running) {
        SHORT start = GetAsyncKeyState(VK_F9);

        if (start) {
            std::cout << "starting to generate data! :-)" << std::endl;
            this->is_running = true;
        }
    }
#endif
    if (!this->is_running)
        return;

    // Check if we want to stop
#if defined(_WIN32) && !defined(NDEBUG)
    SHORT was_pressed = GetAsyncKeyState(VK_F10);

    if (was_pressed) {
        std::cout << "Quit button pressed, stopping data generation " << std::endl;

        std::cout << "Info over the generated photons: " << std::endl;

        int i = 0;
        std::cout << "Channel " << i << ": " << this->n_events[i] << " photons. First: " << this->first_times[i] <<
                     " last: " << this->last_times[i] << std::endl;
        i = 1;
        std::cout << "Channel " << i << ": " << this->n_events[i] << " photons. First: " << this->first_times[i] <<
                     " last: " << this->last_times[i] << std::endl;
        i = 2;
        std::cout << "Channel " << i << ": " << this->n_events[i] << " photons. First: " << this->first_times[i] <<
                     " last: " << this->last_times[i] << std::endl;
        this->is_running = false;
        return;
    }
#endif
    //

    double det_prob = this->detection_probability;

    int64_t p_t = this->last_pulse_time;
    // sync_divider: number of pulses to not write out, + 1
    for (int64_t pulse_t = p_t + this->pulse_period_timeunits;
         pulse_t < p_t + duration;
         pulse_t += this->pulse_period_timeunits) {
        this->sync_counter++;

        double lifetime = 1.0/this->decayrate_inv_seconds;

        this->last_largest_time = pulse_t;
        this->last_pulse_time = pulse_t;

        if (this->sync_counter + 1 <= this->sync_divider) {
            // skip
        } else {
            // Save the pulse event
            if (this->chan0_enabled) {
                timestamps->push_back(pulse_t);
                channel_IDs->push_back(0);

                if (this->n_events.at(0) == 0) {
                    this->first_times[0] = pulse_t;
                }

                this->n_events[0]++;
                this->last_times[0] = pulse_t;
            }
            this->sync_counter = 0;
            last_pulse_published = pulse_t;
        }

        for (int64_t i = 0; i < this->n_emitters; i++) {
            double p_det = gsl_ran_flat(this->r, 0.0, 1.0);
            double dt1 = gsl_ran_exponential(this->r, lifetime); // in seconds
            int64_t dt1_ = (int64_t)round(dt1 / this->timeunit_seconds);

            if (p_det > det_prob) {
                continue;
            }

            if (dt1_ > pulse_t + this->pulse_period_timeunits)
                continue; // Just discard 'pile up' photons

            // Randomly decide the detector
            double detector = gsl_ran_flat(this->r, 0.0, 1.0);

            // TODO: implement detector dead time?
            // Alternate, so we can observe anti-bunching
            if (detector < 0.5 && this->chan1_enabled) {
                channel_IDs->push_back(1);
                timestamps->push_back(pulse_t + dt1_);// + this->chan1_induced_delay);

                if (this->n_events[1] == 0) {
                    this->first_times[1] = pulse_t + dt1_;
                }

                this->n_events[1]++;
                this->last_times[1] = pulse_t + dt1_;

            } else if (detector >= 0.5 && this->chan2_enabled) {
                channel_IDs->push_back(2);
                timestamps->push_back(pulse_t + dt1_);//  + this->chan2_induced_delay);

                if (this->n_events[2] == 0) {
                    this->first_times[2] = pulse_t + dt1_;
                }

                this->n_events[2]++;
                this->last_times[2] = pulse_t + dt1_;
            }
        }
    }
}

uint64_t FakeCommunicator::ReceiveData(std::vector<int64_t>* timestamps,
                         std::vector<uint8_t>* channel_IDs)
{

    if (first_run) {
        first_run = false;
        this->prev_time = std::chrono::high_resolution_clock::now();
        return 0;
    }

    auto curr_time = std::chrono::high_resolution_clock::now();
    double dt_seconds = std::chrono::duration<double, std::milli>(curr_time - this->prev_time).count() / 1000.0;
    dt_seconds *= this->sim_speedup_factor;
    this->prev_time = curr_time;

    int64_t dt_timeunits = 0;

    {
        std::lock_guard<std::mutex> l(this->m);
        dt_timeunits = dt_seconds / this->timeunit_seconds;
    }

    this->gen_pulses(timestamps, channel_IDs, dt_timeunits);
    return timestamps->size();
}
