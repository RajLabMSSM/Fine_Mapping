'use strict';

/**
 * Scatter Data Layer
 * Implements a standard scatter plot
 * @class LocusZoom.DataLayers.scatter
 */
LocusZoom.DataLayers.add('scatter', function(layout) {
    // Define a default layout for this DataLayer type and merge it with the passed argument
    this.DefaultLayout = {
        point_size: 40,
        point_shape: 'circle',
        tooltip_positioning: 'horizontal',
        color: '#888888',
        fill_opacity: 1,
        y_axis: {
            axis: 1
        },
        id_field: 'id'
    };
    layout = LocusZoom.Layouts.merge(layout, this.DefaultLayout);

    // Extra default for layout spacing
    // Not in default layout since that would make the label attribute always present
    if (layout.label && isNaN(layout.label.spacing)) {
        layout.label.spacing = 4;
    }

    // Apply the arguments to set LocusZoom.DataLayer as the prototype
    LocusZoom.DataLayer.apply(this, arguments);

    // Reimplement the positionTooltip() method to be scatter-specific
    this.positionTooltip = function(id) {
        if (typeof id != 'string') {
            throw new Error('Unable to position tooltip: id is not a string');
        }
        if (!this.tooltips[id]) {
            throw new Error('Unable to position tooltip: id does not point to a valid tooltip');
        }
        var top, left, arrow_type, arrow_top, arrow_left;
        var tooltip = this.tooltips[id];
        var point_size = this.resolveScalableParameter(this.layout.point_size, tooltip.data);
        var offset = Math.sqrt(point_size / Math.PI);
        var arrow_width = 7; // as defined in the default stylesheet
        var stroke_width = 1; // as defined in the default stylesheet
        var border_radius = 6; // as defined in the default stylesheet
        var page_origin = this.getPageOrigin();
        var x_center = this.parent.x_scale(tooltip.data[this.layout.x_axis.field]);
        var y_scale  = 'y' + this.layout.y_axis.axis + '_scale';
        var y_center = this.parent[y_scale](tooltip.data[this.layout.y_axis.field]);
        var tooltip_box = tooltip.selector.node().getBoundingClientRect();
        var data_layer_height = this.parent.layout.height - (this.parent.layout.margin.top + this.parent.layout.margin.bottom);
        var data_layer_width = this.parent.layout.width - (this.parent.layout.margin.left + this.parent.layout.margin.right);
        if (this.layout.tooltip_positioning === 'vertical') {
            // Position horizontally centered above the point
            var offset_right = Math.max((tooltip_box.width / 2) - x_center, 0);
            var offset_left = Math.max((tooltip_box.width / 2) + x_center - data_layer_width, 0);
            left = page_origin.x + x_center - (tooltip_box.width / 2) - offset_left + offset_right;
            arrow_left = (tooltip_box.width / 2) - (arrow_width / 2) + offset_left - offset_right - offset;
            // Position vertically above the point unless there's insufficient space, then go below
            if (tooltip_box.height + stroke_width + arrow_width > data_layer_height - (y_center + offset)) {
                top = page_origin.y + y_center - (offset + tooltip_box.height + stroke_width + arrow_width);
                arrow_type = 'down';
                arrow_top = tooltip_box.height - stroke_width;
            } else {
                top = page_origin.y + y_center + offset + stroke_width + arrow_width;
                arrow_type = 'up';
                arrow_top = 0 - stroke_width - arrow_width;
            }
        } else {
            // Position horizontally on the left or the right depending on which side of the plot the point is on
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
            data_layer_height = this.parent.layout.height - (this.parent.layout.margin.top + this.parent.layout.margin.bottom);
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

    // Function to flip labels from being anchored at the start of the text to the end
    // Both to keep labels from running outside the data layer and  also as a first
    // pass on recursive separation
    this.flip_labels = function() {
        var data_layer = this;
        var point_size = data_layer.resolveScalableParameter(data_layer.layout.point_size, {});
        var spacing = data_layer.layout.label.spacing;
        var handle_lines = Boolean(data_layer.layout.label.lines);
        var min_x = 2 * spacing;
        var max_x = data_layer.parent.layout.width - data_layer.parent.layout.margin.left - data_layer.parent.layout.margin.right - (2 * spacing);
        var flip = function(dn, dnl) {
            var dnx = +dn.attr('x');
            var text_swing = (2 * spacing) + (2 * Math.sqrt(point_size));
            if (handle_lines) {
                var dnlx2 = +dnl.attr('x2');
                var line_swing = spacing + (2 * Math.sqrt(point_size));
            }
            if (dn.style('text-anchor') === 'start') {
                dn.style('text-anchor', 'end');
                dn.attr('x', dnx - text_swing);
                if (handle_lines) { dnl.attr('x2', dnlx2 - line_swing); }
            } else {
                dn.style('text-anchor', 'start');
                dn.attr('x', dnx + text_swing);
                if (handle_lines) { dnl.attr('x2', dnlx2 + line_swing); }
            }
        };
        // Flip any going over the right edge from the right side to the left side
        // (all labels start on the right side)
        data_layer.label_texts.each(function (d, i) {
            var a = this;
            var da = d3.select(a);
            var dax = +da.attr('x');
            var abound = da.node().getBoundingClientRect();
            if (dax + abound.width + spacing > max_x) {
                var dal = handle_lines ? d3.select(data_layer.label_lines[0][i]) : null;
                flip(da, dal);
            }
        });
        // Second pass to flip any others that haven't flipped yet if they collide with another label
        data_layer.label_texts.each(function (d, i) {
            var a = this;
            var da = d3.select(a);
            if (da.style('text-anchor') === 'end') {
                return;
            }
            var dax = +da.attr('x');
            var abound = da.node().getBoundingClientRect();
            var dal = handle_lines ? d3.select(data_layer.label_lines[0][i]) : null;
            data_layer.label_texts.each(function () {
                var b = this;
                var db = d3.select(b);
                var bbound = db.node().getBoundingClientRect();
                var collision = abound.left < bbound.left + bbound.width + (2 * spacing) &&
                    abound.left + abound.width + (2 * spacing) > bbound.left &&
                    abound.top < bbound.top + bbound.height + (2 * spacing) &&
                    abound.height + abound.top + (2 * spacing) > bbound.top;
                if (collision) {
                    flip(da, dal);
                    // Double check that this flip didn't push the label past min_x. If it did, immediately flip back.
                    dax = +da.attr('x');
                    if (dax - abound.width - spacing < min_x) {
                        flip(da, dal);
                    }
                }
                return;
            });
        });
    };

    // Recursive function to space labels apart immediately after initial render
    // Adapted from thudfactor's fiddle here: https://jsfiddle.net/thudfactor/HdwTH/
    // TODO: Make labels also aware of data elements
    this.separate_labels = function() {
        this.seperate_iterations++;
        var data_layer = this;
        var alpha = 0.5;
        if (!this.layout.label) {
            // Guard against layout changing in the midst of iterative rerender
            return;
        }
        var spacing = this.layout.label.spacing;
        var again = false;
        data_layer.label_texts.each(function () {
            var a = this;
            var da = d3.select(a);
            var y1 = da.attr('y');
            data_layer.label_texts.each(function () {
                var b = this;
                // a & b are the same element and don't collide.
                if (a === b) {
                    return;
                }
                var db = d3.select(b);
                // a & b are on opposite sides of the chart and
                // don't collide
                if (da.attr('text-anchor') !== db.attr('text-anchor')) {
                    return;
                }
                // Determine if the  bounding rects for the two text elements collide
                var abound = da.node().getBoundingClientRect();
                var bbound = db.node().getBoundingClientRect();
                var collision = abound.left < bbound.left + bbound.width + (2 * spacing) &&
                    abound.left + abound.width + (2 * spacing) > bbound.left &&
                    abound.top < bbound.top + bbound.height + (2 * spacing) &&
                    abound.height + abound.top + (2 * spacing) > bbound.top;
                if (!collision) {
                    return;
                }
                again = true;
                // If the labels collide, we'll push each
                // of the two labels up and down a little bit.
                var y2 = db.attr('y');
                var sign = abound.top < bbound.top ? 1 : -1;
                var adjust = sign * alpha;
                var new_a_y = +y1 - adjust;
                var new_b_y = +y2 + adjust;
                // Keep new values from extending outside the data layer
                var min_y = 2 * spacing;
                var max_y = data_layer.parent.layout.height - data_layer.parent.layout.margin.top - data_layer.parent.layout.margin.bottom - (2 * spacing);
                var delta;
                if (new_a_y - (abound.height / 2) < min_y) {
                    delta = +y1 - new_a_y;
                    new_a_y = +y1;
                    new_b_y += delta;
                } else if (new_b_y - (bbound.height / 2) < min_y) {
                    delta = +y2 - new_b_y;
                    new_b_y = +y2;
                    new_a_y += delta;
                }
                if (new_a_y + (abound.height / 2) > max_y) {
                    delta = new_a_y - +y1;
                    new_a_y = +y1;
                    new_b_y -= delta;
                } else if (new_b_y + (bbound.height / 2) > max_y) {
                    delta = new_b_y - +y2;
                    new_b_y = +y2;
                    new_a_y -= delta;
                }
                da.attr('y',new_a_y);
                db.attr('y',new_b_y);
            });
        });
        if (again) {
            // Adjust lines to follow the labels
            if (data_layer.layout.label.lines) {
                var label_elements = data_layer.label_texts[0];
                data_layer.label_lines.attr('y2',function(d,i) {
                    var label_line = d3.select(label_elements[i]);
                    return label_line.attr('y');
                });
            }
            // After ~150 iterations we're probably beyond diminising returns, so stop recursing
            if (this.seperate_iterations < 150) {
                setTimeout(function() {
                    this.separate_labels();
                }.bind(this), 1);
            }
        }
    };

    // Implement the main render function
    this.render = function() {

        var data_layer = this;
        var x_scale = 'x_scale';
        var y_scale = 'y' + this.layout.y_axis.axis + '_scale';

        if (this.layout.label) {
            // Apply filters to generate a filtered data set
            var filtered_data = this.data.filter(function(d) {
                if (!data_layer.layout.label.filters) {
                    return true;
                } else {
                    // Start by assuming a match, run through all filters to test if not a match on any one
                    var match = true;
                    data_layer.layout.label.filters.forEach(function(filter) {
                        var field_value = (new LocusZoom.Data.Field(filter.field)).resolve(d);
                        if (['!=', '='].indexOf(filter.operator) === -1 && isNaN(field_value)) {
                            // If the filter can only be used with numbers, then the value must be numeric.
                            match = false;
                        } else {
                            switch (filter.operator) {
                            case '<':
                                if (!(field_value < filter.value)) { match = false; }
                                break;
                            case '<=':
                                if (!(field_value <= filter.value)) { match = false; }
                                break;
                            case '>':
                                if (!(field_value > filter.value)) { match = false; }
                                break;
                            case '>=':
                                if (!(field_value >= filter.value)) { match = false; }
                                break;
                            case '=':
                                if (!(field_value === filter.value)) { match = false; }
                                break;
                            case '!=':
                                // Deliberately allow weak comparisons to test for "anything with a value present" (null or undefined)
                                // eslint-disable-next-line eqeqeq
                                if (field_value == filter.value) { match = false; }
                                break;
                            default:
                                // If we got here the operator is not valid, so the filter should fail
                                match = false;
                                break;
                            }
                        }
                    });
                    return match;
                }
            });
            // Render label groups
            var self = this;
            this.label_groups = this.svg.group
                .selectAll('g.lz-data_layer-' + this.layout.type + '-label')
                .data(filtered_data, function(d) { return d[self.layout.id_field]  + '_label'; });
            this.label_groups.enter()
                .append('g')
                .attr('class', 'lz-data_layer-' + this.layout.type + '-label');
            // Render label texts
            if (this.label_texts) { this.label_texts.remove(); }
            this.label_texts = this.label_groups.append('text')
                .attr('class', 'lz-data_layer-' + this.layout.type + '-label');
            this.label_texts
                .text(function(d) {
                    return LocusZoom.parseFields(d, data_layer.layout.label.text || '');
                })
                .style(data_layer.layout.label.style || {})
                .attr({
                    'x': function(d) {
                        var x = data_layer.parent[x_scale](d[data_layer.layout.x_axis.field])
                              + Math.sqrt(data_layer.resolveScalableParameter(data_layer.layout.point_size, d))
                              + data_layer.layout.label.spacing;
                        if (isNaN(x)) { x = -1000; }
                        return x;
                    },
                    'y': function(d) {
                        var y = data_layer.parent[y_scale](d[data_layer.layout.y_axis.field]);
                        if (isNaN(y)) { y = -1000; }
                        return y;
                    },
                    'text-anchor': function() {
                        return 'start';
                    }
                });
            // Render label lines
            if (data_layer.layout.label.lines) {
                if (this.label_lines) { this.label_lines.remove(); }
                this.label_lines = this.label_groups.append('line')
                    .attr('class', 'lz-data_layer-' + this.layout.type + '-label');
                this.label_lines
                    .style(data_layer.layout.label.lines.style || {})
                    .attr({
                        'x1': function(d) {
                            var x = data_layer.parent[x_scale](d[data_layer.layout.x_axis.field]);
                            if (isNaN(x)) { x = -1000; }
                            return x;
                        },
                        'y1': function(d) {
                            var y = data_layer.parent[y_scale](d[data_layer.layout.y_axis.field]);
                            if (isNaN(y)) { y = -1000; }
                            return y;
                        },
                        'x2': function(d) {
                            var x = data_layer.parent[x_scale](d[data_layer.layout.x_axis.field])
                                  + Math.sqrt(data_layer.resolveScalableParameter(data_layer.layout.point_size, d))
                                  + (data_layer.layout.label.spacing / 2);
                            if (isNaN(x)) { x = -1000; }
                            return x;
                        },
                        'y2': function(d) {
                            var y = data_layer.parent[y_scale](d[data_layer.layout.y_axis.field]);
                            if (isNaN(y)) { y = -1000; }
                            return y;
                        }
                    });
            }
            // Remove labels when they're no longer in the filtered data set
            this.label_groups.exit().remove();
        } else {
            // If the layout definition has changed (& no longer specifies labels), strip any previously rendered
            if (this.label_groups) { this.label_groups.remove(); }
            if (this.label_lines) { this.label_lines.remove(); }
        }

        // Generate main scatter data elements
        var selection = this.svg.group
            .selectAll('path.lz-data_layer-' + this.layout.type)
            .data(this.data, function(d) { return d[this.layout.id_field]; }.bind(this));

        // Create elements, apply class, ID, and initial position
        var initial_y = isNaN(this.parent.layout.height) ? 0 : this.parent.layout.height;
        selection.enter()
            .append('path')
            .attr('class', 'lz-data_layer-' + this.layout.type)
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
            selection
                .transition()
                .duration(this.layout.transition.duration || 0)
                .ease(this.layout.transition.ease || 'cubic-in-out')
                .attr('transform', transform)
                .attr('fill', fill)
                .attr('fill-opacity', fill_opacity)
                .attr('d', shape);
        } else {
            selection
                .attr('transform', transform)
                .attr('fill', fill)
                .attr('fill-opacity', fill_opacity)
                .attr('d', shape);
        }

        // Remove old elements as needed
        selection.exit().remove();

        // Apply default event emitters to selection
        selection.on('click.event_emitter', function(element) {
            this.parent.emit('element_clicked', element, true);
        }.bind(this));

        // Apply mouse behaviors
        this.applyBehaviors(selection);

        // Apply method to keep labels from overlapping each other
        if (this.layout.label) {
            this.flip_labels();
            this.seperate_iterations = 0;
            this.separate_labels();
            // Apply default event emitters to selection
            this.label_texts.on('click.event_emitter', function(element) {
                this.parent.emit('element_clicked', element, true);
            }.bind(this));
            // Extend mouse behaviors to labels
            this.applyBehaviors(this.label_texts);
        }

    };

    // Method to set a passed element as the LD reference in the plot-level state
    this.makeLDReference = function(element) {
        var ref = null;
        if (typeof element == 'undefined') {
            throw new Error('makeLDReference requires one argument of any type');
        } else if (typeof element == 'object') {
            if (this.layout.id_field && typeof element[this.layout.id_field] != 'undefined') {
                ref = element[this.layout.id_field].toString();
            } else if (typeof element['id'] != 'undefined') {
                ref = element['id'].toString();
            } else {
                ref = element.toString();
            }
        } else {
            ref = element.toString();
        }
        this.parent_plot.applyState({ ldrefvar: ref });
    };

    return this;

});

/**
 * A scatter plot in which the x-axis represents categories, rather than individual positions.
 * For example, this can be used by PheWAS plots to show related groups. This plot allows the categories to be
 *   determined dynamically when data is first loaded.
 *
 * @class LocusZoom.DataLayers.category_scatter
 * @augments LocusZoom.DataLayers.scatter
 */
LocusZoom.DataLayers.extend('scatter', 'category_scatter', {
    /**
     * This plot layer makes certain assumptions about the data passed in. Transform the raw array of records from
     *   the datasource to prepare it for plotting, as follows:
     * 1. The scatter plot assumes that all records are given in sequence (pre-grouped by `category_field`)
     * 2. It assumes that all records have an x coordinate for individual plotting
     * @private
     */
    _prepareData: function() {
        var xField = this.layout.x_axis.field || 'x';
        // The (namespaced) field from `this.data` that will be used to assign datapoints to a given category & color
        var category_field = this.layout.x_axis.category_field;
        if (!category_field) {
            throw new Error('Layout for ' + this.layout.id + ' must specify category_field');
        }
        // Sort the data so that things in the same category are adjacent (case-insensitive by specified field)
        var sourceData = this.data
            .sort(function(a, b) {
                var ak = a[category_field];
                var bk = b[category_field];
                var av = ak.toString ? ak.toString().toLowerCase() : ak;
                var bv = bk.toString ? bk.toString().toLowerCase() : bk;
                return (av === bv) ? 0 : (av < bv ? -1 : 1);});
        sourceData.forEach(function(d, i) {
            // Implementation detail: Scatter plot requires specifying an x-axis value, and most datasources do not
            //   specify plotting positions. If a point is missing this field, fill in a synthetic value.
            d[xField] = d[xField] || i;
        });
        return sourceData;
    },

    /**
     * Identify the unique categories on the plot, and update the layout with an appropriate color scheme.
     * Also identify the min and max x value associated with the category, which will be used to generate ticks
     * @private
     * @returns {Object.<String, Number[]>} Series of entries used to build category name ticks {category_name: [min_x, max_x]}
     */
    _generateCategoryBounds: function() {
        // TODO: API may return null values in category_field; should we add placeholder category label?
        // The (namespaced) field from `this.data` that will be used to assign datapoints to a given category & color
        var category_field = this.layout.x_axis.category_field;
        var xField = this.layout.x_axis.field || 'x';
        var uniqueCategories = {};
        this.data.forEach(function(item) {
            var category = item[category_field];
            var x = item[xField];
            var bounds = uniqueCategories[category] || [x, x];
            uniqueCategories[category] = [Math.min(bounds[0], x), Math.max(bounds[1], x)];
        });

        var categoryNames = Object.keys(uniqueCategories);
        this._setDynamicColorScheme(categoryNames);

        return uniqueCategories;
    },

    /**
     * Automatically define a color scheme for the layer based on data returned from the server.
     *   If part of the color scheme has been specified, it will fill in remaining missing information.
     *
     * There are three scenarios:
     * 1. The layout does not specify either category names or (color) values. Dynamically build both based on
     *    the data and update the layout.
     * 2. The layout specifies colors, but not categories. Use that exact color information provided, and dynamically
     *     determine what categories are present in the data. (cycle through the available colors, reusing if there
     *     are a lot of categories)
     * 3. The layout specifies exactly what colors and categories to use (and they match the data!). This is useful to
     *    specify an explicit mapping between color scheme and category names, when you want to be sure that the
     *    plot matches a standard color scheme.
     *    (If the layout specifies categories that do not match the data, the user specified categories will be ignored)
     *
     * This method will only act if the layout defines a `categorical_bin` scale function for coloring. It may be
     *   overridden in a subclass to suit other types of coloring methods.
     *
     * @param {String[]} categoryNames
     * @private
     */
    _setDynamicColorScheme: function(categoryNames) {
        var colorParams = this.layout.color.parameters;
        var baseParams = this._base_layout.color.parameters;

        // If the layout does not use a supported coloring scheme, or is already complete, this method should do nothing
        if (this.layout.color.scale_function !== 'categorical_bin') {
            throw new Error('This layer requires that coloring be specified as a `categorical_bin`');
        }

        if (baseParams.categories.length && baseParams.values.length) {
            // If there are preset category/color combos, make sure that they apply to the actual dataset
            var parameters_categories_hash = {};
            baseParams.categories.forEach(function (category) { parameters_categories_hash[category] = 1; });
            if (categoryNames.every(function (name) { return parameters_categories_hash.hasOwnProperty(name); })) {
                // The layout doesn't have to specify categories in order, but make sure they are all there
                colorParams.categories = baseParams.categories;
            } else {
                colorParams.categories = categoryNames;
            }
        } else {
            colorParams.categories = categoryNames;
        }
        // Prefer user-specified colors if provided. Make sure that there are enough colors for all the categories.
        var colors;
        if (baseParams.values.length) {
            colors = baseParams.values;
        } else {
            var color_scale = categoryNames.length <= 10 ? d3.scale.category10 : d3.scale.category20;
            colors = color_scale().range();
        }
        while (colors.length < categoryNames.length) { colors = colors.concat(colors); }
        colors = colors.slice(0, categoryNames.length);  // List of hex values, should be of same length as categories array
        colorParams.values = colors;
    },

    /**
     *
     * @param dimension
     * @param {Object} [config] Parameters that customize how ticks are calculated (not style)
     * @param {('left'|'center'|'right')} [config.position='left'] Align ticks with the center or edge of category
     * @returns {Array}
     */
    getTicks: function(dimension, config) { // Overrides parent method
        if (['x', 'y1', 'y2'].indexOf(dimension) === -1) {
            throw new Error('Invalid dimension identifier');
        }
        var position = config.position || 'left';
        if (['left', 'center', 'right'].indexOf(position) === -1) {
            throw new Error('Invalid tick position');
        }

        var categoryBounds = this._categories;
        if (!categoryBounds || !Object.keys(categoryBounds).length) {
            return [];
        }

        if (dimension === 'y') {
            return [];
        }

        if (dimension === 'x') {
            // If colors have been defined by this layer, use them to make tick colors match scatterplot point colors
            var knownCategories = this.layout.color.parameters.categories || [];
            var knownColors = this.layout.color.parameters.values || [];

            return Object.keys(categoryBounds).map(function (category, index) {
                var bounds = categoryBounds[category];
                var xPos;

                switch(position) {
                case 'left':
                    xPos = bounds[0];
                    break;
                case 'center':
                    // Center tick under one or many elements as appropriate
                    var diff = bounds[1] - bounds[0];
                    xPos = bounds[0] + (diff !== 0 ? diff : bounds[0]) / 2;
                    break;
                case 'right':
                    xPos = bounds[1];
                    break;
                }
                return {
                    x: xPos,
                    text: category,
                    style: {
                        'fill': knownColors[knownCategories.indexOf(category)] || '#000000'
                    }
                };
            });
        }
    },

    applyCustomDataMethods: function() {
        this.data = this._prepareData();
        /**
         * Define category names and extents (boundaries) for plotting.  TODO: properties in constructor
         * @member {Object.<String, Number[]>} Category names and extents, in the form {category_name: [min_x, max_x]}
         */
        this._categories = this._generateCategoryBounds();
        return this;
    }
});
