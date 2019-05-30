'use strict';

/**
 * Forest Data Layer
 * Implements a standard forest plot. In order to space out points, any layout using this must specify axis ticks
 *  and extent in advance.
 *
 * If you are using dynamically fetched data, consider using `category_forest` instead.
 *
 * @class LocusZoom.DataLayers.forest
 */
LocusZoom.DataLayers.add('forest', function(layout) {

    // Define a default layout for this DataLayer type and merge it with the passed argument
    this.DefaultLayout = {
        point_size: 40,
        point_shape: 'square',
        color: '#888888',
        fill_opacity: 1,
        y_axis: {
            axis: 2
        },
        id_field: 'id',
        confidence_intervals: {
            start_field: 'ci_start',
            end_field: 'ci_end'
        },
        show_no_significance_line: true
    };
    layout = LocusZoom.Layouts.merge(layout, this.DefaultLayout);

    // Apply the arguments to set LocusZoom.DataLayer as the prototype
    LocusZoom.DataLayer.apply(this, arguments);

    // Reimplement the positionTooltip() method to be forest-specific
    this.positionTooltip = function(id) {
        if (typeof id != 'string') {
            throw new Error('Unable to position tooltip: id is not a string');
        }
        if (!this.tooltips[id]) {
            throw new Error('Unable to position tooltip: id does not point to a valid tooltip');
        }
        var tooltip = this.tooltips[id];
        var point_size = this.resolveScalableParameter(this.layout.point_size, tooltip.data);
        var arrow_width = 7; // as defined in the default stylesheet
        var stroke_width = 1; // as defined in the default stylesheet
        var border_radius = 6; // as defined in the default stylesheet
        var page_origin = this.getPageOrigin();
        var x_center = this.parent.x_scale(tooltip.data[this.layout.x_axis.field]);
        var y_scale  = 'y' + this.layout.y_axis.axis + '_scale';
        var y_center = this.parent[y_scale](tooltip.data[this.layout.y_axis.field]);
        var tooltip_box = tooltip.selector.node().getBoundingClientRect();
        // Position horizontally on the left or the right depending on which side of the plot the point is on
        var offset = Math.sqrt(point_size / Math.PI);
        var left, arrow_type, arrow_left;
        if (x_center <= this.parent.layout.width / 2) {
            left = page_origin.x + x_center + offset + arrow_width + stroke_width;
            arrow_type = 'left';
            arrow_left = -1 * (arrow_width + stroke_width);
        } else {
            left = page_origin.x + x_center - tooltip_box.width - offset - arrow_width - stroke_width;
            arrow_type = 'right';
            arrow_left = tooltip_box.width - stroke_width;
        }
        // Position vertically centered unless we're at the top or bottom of the plot
        var data_layer_height = this.parent.layout.height - (this.parent.layout.margin.top + this.parent.layout.margin.bottom);
        var top, arrow_top;
        if (y_center - (tooltip_box.height / 2) <= 0) { // Too close to the top, push it down
            top = page_origin.y + y_center - (1.5 * arrow_width) - border_radius;
            arrow_top = border_radius;
        } else if (y_center + (tooltip_box.height / 2) >= data_layer_height) { // Too close to the bottom, pull it up
            top = page_origin.y + y_center + arrow_width + border_radius - tooltip_box.height;
            arrow_top = tooltip_box.height - (2 * arrow_width) - border_radius;
        } else { // vertically centered
            top = page_origin.y + y_center - (tooltip_box.height / 2);
            arrow_top = (tooltip_box.height / 2) - arrow_width;
        }
        // Apply positions to the main div
        tooltip.selector.style('left', left + 'px').style('top', top + 'px');
        // Create / update position on arrow connecting tooltip to data
        if (!tooltip.arrow) {
            tooltip.arrow = tooltip.selector.append('div').style('position', 'absolute');
        }
        tooltip.arrow
            .attr('class', 'lz-data_layer-tooltip-arrow_' + arrow_type)
            .style('left', arrow_left + 'px')
            .style('top', arrow_top + 'px');
    };

    // Implement the main render function
    this.render = function() {

        var x_scale = 'x_scale';
        var y_scale = 'y' + this.layout.y_axis.axis + '_scale';

        // Generate confidence interval paths if fields are defined
        if (this.layout.confidence_intervals
            && this.layout.fields.indexOf(this.layout.confidence_intervals.start_field) !== -1
            && this.layout.fields.indexOf(this.layout.confidence_intervals.end_field) !== -1) {
            // Generate a selection for all forest plot confidence intervals
            var ci_selection = this.svg.group
                .selectAll('rect.lz-data_layer-forest.lz-data_layer-forest-ci')
                .data(this.data, function(d) { return d[this.layout.id_field]; }.bind(this));
            // Create confidence interval rect elements
            ci_selection.enter()
                .append('rect')
                .attr('class', 'lz-data_layer-forest lz-data_layer-forest-ci')
                .attr('id', function(d) { return this.getElementId(d) + '_ci'; }.bind(this))
                .attr('transform', 'translate(0,' + (isNaN(this.parent.layout.height) ? 0 : this.parent.layout.height) + ')');
            // Apply position and size parameters using transition if necessary
            var ci_transform = function(d) {
                var x = this.parent[x_scale](d[this.layout.confidence_intervals.start_field]);
                var y = this.parent[y_scale](d[this.layout.y_axis.field]);
                if (isNaN(x)) { x = -1000; }
                if (isNaN(y)) { y = -1000; }
                return 'translate(' + x + ',' + y + ')';
            }.bind(this);
            var ci_width = function(d) {
                return this.parent[x_scale](d[this.layout.confidence_intervals.end_field])
                     - this.parent[x_scale](d[this.layout.confidence_intervals.start_field]);
            }.bind(this);
            var ci_height = 1;
            if (this.canTransition()) {
                ci_selection
                    .transition()
                    .duration(this.layout.transition.duration || 0)
                    .ease(this.layout.transition.ease || 'cubic-in-out')
                    .attr('transform', ci_transform)
                    .attr('width', ci_width).attr('height', ci_height);
            } else {
                ci_selection
                    .attr('transform', ci_transform)
                    .attr('width', ci_width).attr('height', ci_height);
            }
            // Remove old elements as needed
            ci_selection.exit().remove();
        }

        // Generate a selection for all forest plot points
        var points_selection = this.svg.group
            .selectAll('path.lz-data_layer-forest.lz-data_layer-forest-point')
            .data(this.data, function(d) { return d[this.layout.id_field]; }.bind(this));

        // Create elements, apply class, ID, and initial position
        var initial_y = isNaN(this.parent.layout.height) ? 0 : this.parent.layout.height;
        points_selection.enter()
            .append('path')
            .attr('class', 'lz-data_layer-forest lz-data_layer-forest-point')
            .attr('id', function(d) { return this.getElementId(d); }.bind(this))
            .attr('transform', 'translate(0,' + initial_y + ')');

        // Generate new values (or functions for them) for position, color, size, and shape
        var transform = function(d) {
            var x = this.parent[x_scale](d[this.layout.x_axis.field]);
            var y = this.parent[y_scale](d[this.layout.y_axis.field]);
            if (isNaN(x)) { x = -1000; }
            if (isNaN(y)) { y = -1000; }
            return 'translate(' + x + ',' + y + ')';
        }.bind(this);

        var fill = function(d) { return this.resolveScalableParameter(this.layout.color, d); }.bind(this);
        var fill_opacity = function(d) { return this.resolveScalableParameter(this.layout.fill_opacity, d); }.bind(this);

        var shape = d3.svg.symbol()
            .size(function(d) { return this.resolveScalableParameter(this.layout.point_size, d); }.bind(this))
            .type(function(d) { return this.resolveScalableParameter(this.layout.point_shape, d); }.bind(this));

        // Apply position and color, using a transition if necessary
        if (this.canTransition()) {
            points_selection
                .transition()
                .duration(this.layout.transition.duration || 0)
                .ease(this.layout.transition.ease || 'cubic-in-out')
                .attr('transform', transform)
                .attr('fill', fill)
                .attr('fill-opacity', fill_opacity)
                .attr('d', shape);
        } else {
            points_selection
                .attr('transform', transform)
                .attr('fill', fill)
                .attr('fill-opacity', fill_opacity)
                .attr('d', shape);
        }

        // Remove old elements as needed
        points_selection.exit().remove();

        // Apply default event emitters to selection
        points_selection.on('click.event_emitter', function(element_data) {
            this.parent.emit('element_clicked', element_data, true);
        }.bind(this));

        // Apply behaviors to points
        this.applyBehaviors(points_selection);
    };
    return this;

});

/**
 * A y-aligned forest plot in which the y-axis represents item labels, which are dynamically chosen when data is loaded.
 *   Each item is assumed to include both data and confidence intervals.
 *   This allows generating forest plots without defining the layout in advance.
 *
 * @class LocusZoom.DataLayers.category_forest
 * @augments LocusZoom.DataLayers.forest
 */
LocusZoom.DataLayers.extend('forest', 'category_forest', {
    _getDataExtent: function(data, axis_config) {
        // In a forest plot, the data range is determined by *three* fields (beta + CI start/end)
        var ci_config = this.layout.confidence_intervals;
        if (ci_config
            && this.layout.fields.indexOf(ci_config.start_field) !== -1
            && this.layout.fields.indexOf(ci_config.end_field) !== -1) {
            var min = function(d) {
                var f = new LocusZoom.Data.Field(ci_config.start_field);
                return +f.resolve(d);
            };

            var max = function(d) {
                var f = new LocusZoom.Data.Field(ci_config.end_field);
                return +f.resolve(d);
            };

            return [d3.min(data, min), d3.max(data, max)];
        }

        // If there are no confidence intervals set, then range must depend only on a single field
        return LocusZoom.DataLayer.prototype._getDataExtent.call(this, data, axis_config);
    },

    getTicks: function(dimension, config) { // Overrides parent method
        if (['x', 'y1', 'y2'].indexOf(dimension) === -1) {
            throw new Error('Invalid dimension identifier' + dimension);
        }

        // Design assumption: one axis (y1 or y2) has the ticks, and the layout says which to use
        // Also assumes that every tick gets assigned a unique matching label
        var axis_num = this.layout.y_axis.axis;
        if (dimension === ('y' + axis_num)) {
            var category_field = this.layout.y_axis.category_field;
            if (!category_field) {
                throw new Error('Layout for ' + this.layout.id + ' must specify category_field');
            }

            return this.data.map(function (item, index) {
                return {
                    y: index + 1,
                    text: item[category_field]
                };
            });
        } else {
            return [];
        }
    },

    applyCustomDataMethods: function () {
        // Add a synthetic yaxis field to ensure data is spread out on plot. Then, set axis floor and ceiling to
        //  correct extents.
        var field_to_add = this.layout.y_axis.field;
        if (!field_to_add) {
            throw new Error('Layout for ' + this.layout.id + ' must specify yaxis.field');
        }

        this.data = this.data.map(function (item, index) {
            item[field_to_add] = index + 1;
            return item;
        });
        // Update axis extents based on one label for every point (with a bit of padding above and below)
        this.layout.y_axis.floor = 0;
        this.layout.y_axis.ceiling = this.data.length + 1;
        return this;
    }
});
