/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#ifndef DECAYCURVE_CALC_H
#define DECAYCURVE_CALC_H

#include <phast_gui/base/base_ui_calc.h>
#include <phast_gui/support/typedefs.h>

class ISettingsContainer;

class decaycurve_calc : public Base_ui_calc
{
private:
    std::vector<int64_t> xs;
    std::vector<int64_t> ys;

    double bin_low;
    double bin_high;

    int64_t max_num_histogram_bins;
protected:
    virtual chan_data_map get_times_to_save(DetectionUpdate event) const;

    virtual DataView::data_update arrange_update(DetectionUpdate event, bool rebin_only);

public:
    decaycurve_calc(const decaycurve_calc&) = delete;
    decaycurve_calc& operator=(const decaycurve_calc&) = delete;
    decaycurve_calc(const std::vector<chan_id>& chans_to_sum,
                    uint64_t bin_width,
                    opaque_ptr our_ref,
                    double bin_low,
                    double bin_high,
                    bool display_only_initial_data,
                    double time_unit,
                    int64_t max_num_histogram_bins=1000000);

    virtual void ClearData();
};

#endif // DECAYCURVE_CALC_H
