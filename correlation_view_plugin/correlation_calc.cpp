/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#include "correlation_calc.h"
#include <phast_gui/support/detectionupdate.h>
#include <libtimetag/algos.h>
#include <math.h>

correlation_calc::correlation_calc(const std::vector<chan_id>& chans_to_correlate,
                                   uint64_t bin_width,
                                   opaque_ptr our_ref,
                                   double bin_low,
                                   double bin_high,
                                   bool display_only_initial_data,
                                   double time_unit,
                                   bool correlate_FCS,
                                   int64_t max_num_bins,
                                   int64_t num_lin_FCS_bins) :
    Base_ui_calc(chans_to_correlate, bin_width, our_ref, display_only_initial_data, time_unit),
    leftchannel(chans_to_correlate.at(0)),
    rightchannel(chans_to_correlate.at(1)),
    t_min(INT64_MAX),
    t_max(0),
    bin_low(bin_low),
    bin_high(bin_high),
    n_photons_left(0),
    n_photons_right(0),
    correlate_FCS(correlate_FCS),
    max_num_bins(max_num_bins),
    lin_FCS(false),
    num_lin_FCS_bins(num_lin_FCS_bins)
{
    if (!correlate_FCS) {
        // Check if our length is feasible at all
        int64_t len = linspace_len(static_cast<timestamp>(this->bin_low),
                                   static_cast<timestamp>(this->bin_high), 1, 1, 1);

        // The length may be so large, that it overflows an int64_t...
        uint64_t temp_len = (uint64_t)(this->bin_high - this->bin_low);

        if (len > this->max_num_bins || temp_len > (INT64_MAX-1)) {
            this->bin_width = 1;
            // For now: switch to "FCS-style" correlation,
            // but using linear bins
            this->correlate_FCS = true;
            this->lin_FCS = true;
            // check the length
            double tot_width = this->bin_high - this->bin_low;
            double est_bin_width = tot_width / this->num_lin_FCS_bins;
            int64_t est_bin_width_ = (int64_t)est_bin_width;
            int64_t new_num_bins = linspace_len(static_cast<timestamp>(this->bin_low),
                                                static_cast<timestamp>(this->bin_high),
                                                est_bin_width_, 1, 1);

            this->xs.resize(new_num_bins,0);
            this->ys.resize(new_num_bins-1,0);

            linspace(this->bin_low, this->bin_high, est_bin_width_, 1, 1, this->xs.data(), this->xs.size());

            this->ClearData();
        } else {

            this->xs.resize(len, 0);
            this->ys.resize(len - 1, 0);

            linspace(this->bin_low, this->bin_high, 1, 1, 1, this->xs.data(), this->xs.size());
        }
    } else {
        this->xs.resize(this->bin_width, 0);
        this->ys.resize(this->bin_width - 1, 0);

        std::vector<double> temp(this->bin_width, 0);

        logspace(log10(this->bin_low),
                 log10(this->bin_high),
                 this->bin_width,
                 10,
                 temp.data());

        for (uint64_t i = 0; i < temp.size(); i++) {
            this->xs[i] = (int64_t)temp[i];
        }
    }
}

void correlation_calc::ClearData()
{
    this->ys.clear();
    this->ys.resize(this->xs.size() - 1, 0);

    this->n_photons_left = 0;
    this->n_photons_right = 0;

    this->t_min = INT64_MAX;
    this->t_max = 0;
    data_displayed = false;
}

chan_data_map correlation_calc::get_times_to_save(DetectionUpdate event) const
{
    int64_t lower_bound = this->bin_low;
    int64_t abs_lower_bound = this->t_max + lower_bound;

    chan_data_map ret;

    for (auto pair : event.updates) { // For each channel...

        chan_id chan = pair.first;
        chan_det_update update = pair.second;

        // If we are not interested in this channel,
        // signal that everything can be trashed.
        if (chan != this->leftchannel && chan != this->rightchannel) {

            ret[chan] = update.macro_end;
            continue;
        }
        uint64_t counter = 0;

        for (data_type it = update.macro_start; it != update.macro_end; ++it) {
            counter++;

            if (*it > abs_lower_bound) {
                ret[chan] = it;
                break;
            }
        }
    }

    return ret;
}

void correlation_calc::update_tmin(DetectionUpdate event)
{
    if (event.updates.count(this->leftchannel) == 0 ||
            event.updates.count(this->rightchannel) == 0) {
        // We don't have any data to work with.
        return;
    }

    chan_det_update left = event.updates.at(this->leftchannel);
    chan_det_update right = event.updates.at(this->rightchannel);

    timestamp t_min_= INT64_MAX;

    if (left.n_new_macro != 0) {
        if (*left.new_macro_start < t_min_) {
            t_min_ = *left.new_macro_start;
        }
    }

    if (right.n_new_macro != 0) {
        if (*right.new_macro_start < t_min_) {
            t_min_ = *right.new_macro_start;
        }
    }

    if (t_min_ < this->t_min) {
        this->t_min = t_min_;
    }
}

void correlation_calc::update_tmax(DetectionUpdate event)
{
    if (event.updates.count(this->leftchannel) == 0 ||
            event.updates.count(this->rightchannel) == 0) {
        // We don't have any data to work with.
        return;
    }

    chan_det_update left = event.updates.at(this->leftchannel);
    chan_det_update right = event.updates.at(this->rightchannel);

    timestamp t_max_= 0;

    if (left.n_new_macro != 0) {
        if (*(left.macro_end-1) > t_max_) {
            t_max_ = *(left.macro_end-1);
        }
    }

    if (right.n_new_macro != 0) {
        if (*(right.macro_end-1) > t_max_) {
            t_max_ = *(right.macro_end-1);
        }
    }

    if (t_max_ > this->t_max) {
        this->t_max = t_max_;
    }
}

void correlation_calc::update_n_photons_leftchan(DetectionUpdate event)
{
    if (event.updates.count(this->leftchannel) == 0)
        return;

    this->n_photons_left += event.updates.at(this->leftchannel).n_new_macro;
}

void correlation_calc::update_n_photons_rightchan(DetectionUpdate event)
{
    if (event.updates.count(this->rightchannel) == 0)
        return;

    this->n_photons_right += event.updates.at(this->rightchannel).n_new_macro;
}

DataView::data_update correlation_calc::arrange_update(DetectionUpdate event, bool rebin_only)
{
    DataView::data_update ret;

    if (event.updates.count(this->leftchannel) == 0 ||
            event.updates.count(this->rightchannel) == 0) {
        // We don't have any data to work with.
        return ret;
    }

    chan_det_update left = event.updates.at(this->leftchannel);
    chan_det_update right = event.updates.at(this->rightchannel);

    if (!rebin_only) {
        // See if we have any photons

        this->update_tmin(event);
        this->update_tmax(event);

        this->update_n_photons_leftchan(event);
        this->update_n_photons_rightchan(event);

        // First we correlate *just* the new data on channel 1 to all the data on chan 2
        //  Correlated are then:
        //  old chan 1 <-> old chan 2
        //  new chan 1 <-> old chan 2
        //  new chan 1 <-> new chan 2
        //
        // Still missing is:
        // old chan 1 <-> new chan 2
        //
        // This last part we do in a second step.

        if (this->correlate_FCS) {
            uint64_t n_left_photons = left.n_new_macro;
            uint64_t n_right_photons = right.n_macro;

            if (n_left_photons != 0 && n_right_photons != 0) {
                correlate_many_per_bin(this->xs.data(),
                                       this->xs.size(),
                                       &*left.new_macro_start,
                                       n_left_photons,
                                       &*right.macro_start,
                                       n_right_photons,
                                       this->ys.data(),
                                       this->ys.size());
            }

            n_left_photons = left.n_old_macro;
            n_right_photons = right.n_new_macro;

            if (n_left_photons != 0 && n_right_photons != 0) {
                correlate_many_per_bin(this->xs.data(),
                                       this->xs.size(),
                                       &*left.macro_start,
                                       n_left_photons,
                                       &*right.new_macro_start,
                                       n_right_photons,
                                       this->ys.data(),
                                       this->ys.size());
            }
        } else {
            correlate_unit_bins(this->xs.data(),
                                this->xs.size(),
                                &*left.new_macro_start,
                                left.n_new_macro,
                                &*right.macro_start,
                                right.n_macro,
                                this->ys.data(),
                                this->ys.size());

            correlate_unit_bins(this->xs.data(),
                                this->xs.size(),
                                &*left.macro_start,
                                left.n_old_macro,
                                &*right.new_macro_start,
                                right.n_new_macro,
                                this->ys.data(),
                                this->ys.size());
        }
    }

    std::vector<int64_t> rebinned_xs(rebin_bin_edges_len(this->xs.size(), this->bin_width), 0);
    std::vector<int64_t> rebinned_ys(rebin_len(this->ys.size(), this->bin_width), 0);

    if (this->bin_width > 1 && (!this->correlate_FCS || this->lin_FCS)) {
        rebin(this->ys.data(),
              this->ys.size(),
              this->bin_width,
              rebinned_ys.data(),
              rebinned_ys.size());

        rebin_bin_edges(this->xs.data(),
                        this->xs.size(),
                        this->bin_width,
                        rebinned_xs.data(),
                        rebinned_xs.size());
    } else {
        rebinned_xs = this->xs;
        rebinned_ys = this->ys;
    }

    std::vector<double> xs2;
    std::vector<double> ys2(rebinned_ys.size(), 0.0);
    normalize_correlation(rebinned_ys.data(),
                        rebinned_ys.size(),
                        rebinned_xs.data(),
                        rebinned_xs.size(),
                        this->t_min,
                        this->t_max,
                        this->n_photons_left,
                        this->n_photons_right,
                          ys2.data());


    for (uint64_t i = 0; i < rebinned_xs.size() - 1; i++) {
        xs2.push_back((double)rebinned_xs[i]*this->time_unit*1e9); // TODO: do not hardcode the timeunit
    }

    ret.xs = xs2;
    ret.ys = ys2;

    return ret;
}
