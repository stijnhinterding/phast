/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#ifndef FAKECOMMUNICATOR_H
#define FAKECOMMUNICATOR_H

#include <phast_gui/interfaces/itimetaggercommunicator.h>

#include <gsl/gsl_rng.h>
#include <chrono>
#include <mutex>

class FakeCommunicator : public ITimeTaggerCommunicator
{
private:
    gsl_rng* r;
    double timeunit_seconds;
    double decayrate_inv_seconds;
    double pulse_period_seconds;
    double detection_probability;

    std::chrono::time_point<std::chrono::high_resolution_clock> prev_time;
    bool first_run;

    int64_t last_largest_time;
    int64_t last_pulse_time;
    int64_t pulse_period_timeunits;
    int64_t last_pulse_published;
    uint64_t sync_counter;
    uint64_t sync_divider;

    int64_t chan1_induced_delay;
    int64_t chan2_induced_delay;

    std::string descriptor;

    bool chan0_enabled;
    bool chan1_enabled;
    bool chan2_enabled;

    bool is_running;

    double sim_speedup_factor;
    std::mutex m;

    std::map<chan_id,timestamp> first_times;
    std::map<chan_id,timestamp> last_times;
    std::map<chan_id,uint64_t> n_events;

    int64_t n_emitters;
private:
    void gen_pulses(std::vector<int64_t> *timestamps,
                    std::vector<uint8_t> *channel_IDs, int64_t duration);

public:
    FakeCommunicator(const FakeCommunicator&) = delete;
    FakeCommunicator& operator=(const FakeCommunicator&) = delete;
    ~FakeCommunicator() override;

    FakeCommunicator(double timeunit_seconds,
                     double decayrate_inv_seconds,
                     double pulse_period_seconds,
                     double detection_probability,
                     uint64_t sync_divider=128,
                     int64_t chan1_induced_delay=0,
                     int64_t chan2_induced_delay=0,
                     int64_t n_emitters=1,
                     double sim_speedup_factor=1);

    void SetSimSpeedupFactor(double value)
    {
        std::lock_guard<std::mutex> l(this->m);
        if (value <= 0.001)
            value = 0.001;

        this->sim_speedup_factor = value;
    }

    double GetSimSpeedupFactor() const
    {
        return this->sim_speedup_factor;
    }

    void SetDecayRate(double value)
    {
        std::lock_guard<std::mutex> l(this->m);
        this->decayrate_inv_seconds = value;
    }

    int64_t GetNumEmitters() const
    {
        return this->n_emitters;
    }

    double GetDecayRate() const
    {
        return this->decayrate_inv_seconds;
    }

    void SetNumEmitters(int64_t value)
    {
        std::lock_guard<std::mutex> l(this->m);

        if (value == 0)
            value = 1;

        this->n_emitters = value;
    }
    void SetPulsePeriod(double value)
    {
        std::lock_guard<std::mutex> l(this->m);
        this->pulse_period_seconds = value;
        this->pulse_period_timeunits = this->pulse_period_seconds/this->timeunit_seconds;
    }

    double GetPulsePeriod() const
    {
        return this->pulse_period_seconds;
    }

    void SetDetectionProbability(double value)
    {
        std::lock_guard<std::mutex> l(this->m);
        if (value > 1.0)
            value = 1.0;

        if (value < 0.0)
            value = 0.0;

        this->detection_probability = value;
    }

    double GetDetectionProbability() const
    {
        return this->detection_probability;
    }

    void SetSyncDivider(uint64_t value)
    {
        std::lock_guard<std::mutex> l(this->m);
        this->sync_divider = value;
    }

    virtual uint64_t ReceiveData(std::vector<int64_t>* timestamps,
                             std::vector<uint8_t>* channel_IDs) override;

    virtual uint64_t GetSyncDivider(chan_id channel_id) override
    {
        if (channel_id == 0)
            return this->sync_divider;

        return 1;
    }

    virtual void DisableChan(uint64_t ID);

    virtual bool ConnectedToDevice() override
    {
        return true;
    }

    virtual const std::string& DeviceDescriptor() const override
    {
        return this->descriptor;
    }

    virtual bool TryToConnect(int64_t) override
    {
        return true;
    }

    virtual bool DataLossSinceLastCall() override
    {
        return false;
    }

    virtual double TimeUnit() const override
    {
        return this->timeunit_seconds;
    }

    virtual bool AnyDeviceAvailable() const override
    {
        return true;
    }

    virtual bool Disconnect() override
    {
        return true;
    }

    virtual uint64_t GetNumDevicesConnected() override
    {
        return 1;
    }

    virtual bool IsRealDevice() const
    {
        return false;
    }
};

#endif // FAKECOMMUNICATOR_H
