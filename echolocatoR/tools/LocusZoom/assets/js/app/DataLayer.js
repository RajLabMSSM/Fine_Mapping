/* global LocusZoom */
'use strict';

/**
 * A data layer is an abstract class representing a data set and its graphical representation within a panel
 * @public
 * @class
 * @param {Object} layout A JSON-serializable object describing the layout for this layer
 * @param {LocusZoom.DataLayer|LocusZoom.Panel} parent Where this layout is used
*/
LocusZoom.DataLayer = function(layout, parent) {
    /** @member {Boolean} */
    this.initialized = false;
    /** @member {Number} */
    this.layout_idx = null;

    /** @member {String} */
    this.id     = null;
    /** @member {LocusZoom.Panel} */
    this.parent = parent || null;
    /**
     * @member {{group: d3.selection, container: d3.selection, clipRect: d3.selection}}
     */
    this.svg    = {};

    /** @member {LocusZoom.Plot} */
    this.parent_plot = null;
    if (typeof parent != 'undefined' && parent instanceof LocusZoom.Panel) { this.parent_plot = parent.parent; }

    /** @member {Object} */
    this.layout = LocusZoom.Layouts.merge(layout || {}, LocusZoom.DataLayer.DefaultLayout);
    if (this.layout.id) { this.id = this.layout.id; }

    // Ensure any axes defined in the layout have an explicit axis number (default: 1)
    if (this.layout.x_axis !== {} && typeof this.layout.x_axis.axis !== 'number') { this.layout.x_axis.axis = 1; }
    if (this.layout.y_axis !== {} && typeof this.layout.y_axis.axis !== 'number') { this.layout.y_axis.axis = 1; }

    /**
     * Values in the layout object may change during rendering etc. Retain a copy of the original data layer state
     * @member {Object}
     */
    this._base_layout = JSON.parse(JSON.stringify(this.layout));

    /** @member {Object} */
    this.state = {};
    /** @member {String} */
    this.state_id = null;

    this.setDefaultState();

    // Initialize parameters for storing data and tool tips
    /** @member {Array} */
    this.data = [];
    if (this.layout.tooltip) {
        /** @member {Object} */
        this.tooltips = {};
    }

    // Initialize flags for tracking global statuses
    this.global_statuses = {
        'highlighted': false,
        'selected': false,
        'faded': false,
        'hidden': false
    };

    return this;

};

/**
 * Instruct this datalayer to begin tracking additional fields from data sources (does not guarantee that such a field actually exists)
 *
 * Custom plots can use this to dynamically extend datalayer functionality after the plot is drawn
 *
 *  (since removing core fields may break layer functionality, there is presently no hook for the inverse behavior)
 * @param fieldName
 * @param namespace
 * @param {String|String[]} transformations The name (or array of names) of transformations to apply to this field
 * @returns {String} The raw string added to the fields array
 */
LocusZoom.DataLayer.prototype.addField = function(fieldName, namespace, transformations) {
    if (!fieldName || !namespace) {
        throw new Error('Must specify field name and namespace to use when adding field');
    }
    var fieldString = namespace + ':' + fieldName;
    if (transformations) {
        fieldString += '|';
        if (typeof transformations === 'string') {
            fieldString += transformations;
        } else if (Array.isArray(transformations)) {
            fieldString += transformations.join('|');
        } else {
            throw new Error('Must provide transformations as either a string or array of strings');
        }
    }
    var fields = this.layout.fields;
    if (fields.indexOf(fieldString) === -1) {
        fields.push(fieldString);
    }
    return fieldString;
};

/**
 * Define default state that should get tracked during the lifetime of this layer.
 *
 * In some special custom usages, it may be useful to completely reset a panel (eg "click for
 *   genome region" links), plotting new data that invalidates any previously tracked state.  This hook makes it
 *   possible to reset without destroying the panel entirely. It is used by `Plot.clearPanelData`.
 */
LocusZoom.DataLayer.prototype.setDefaultState = function() {
    // Define state parameters specific to this data layer. Within plot state, this will live under a key
    //  `panel_name.layer_name`.
    if (this.parent) {
        this.state = this.parent.state;
        this.state_id = this.parent.id + '.' + this.id;
        this.state[this.state_id] = this.state[this.state_id] || {};
        LocusZoom.DataLayer.Statuses.adjectives.forEach(function(status) {
            this.state[this.state_id][status] = this.state[this.state_id][status] || [];
        }.bind(this));
    }
};

/**
 * A basic description of keys expected in a layout. Not intended to be directly used or modified by an end user.
 * @protected
 * @type {{type: string, fields: Array, x_axis: {}, y_axis: {}}}
 */
LocusZoom.DataLayer.DefaultLayout = {
    type: '',
    fields: [],
    x_axis: {},
    y_axis: {}
};

/**
 * Available statuses that individual elements can have. Each status is described by
 *   a verb/antiverb and an adjective. Verbs and antiverbs are used to generate data layer
 *   methods for updating the status on one or more elements. Adjectives are used in class
 *   names and applied or removed from elements to have a visual representation of the status,
 *   as well as used as keys in the state for tracking which elements are in which status(es)
 * @static
 * @type {{verbs: String[], adjectives: String[], menu_antiverbs: String[]}}
 */
LocusZoom.DataLayer.Statuses = {
    verbs: ['highlight', 'select', 'fade', 'hide'],
    adjectives: ['highlighted', 'selected', 'faded', 'hidden'],
    menu_antiverbs: ['unhighlight', 'deselect', 'unfade', 'show']
};

/**
 * Get the fully qualified identifier for the data layer, prefixed by any parent or container elements
 *
 * @returns {string} A dot-delimited string of the format <plot>.<panel>.<data_layer>
 */
LocusZoom.DataLayer.prototype.getBaseId = function() {
    return this.parent_plot.id + '.' + this.parent.id + '.' + this.id;
};

/**
 * Determine the pixel height of data-bound objects represented inside this data layer. (excluding elements such as axes)
 *
 * May be used by operations that resize the data layer to fit available data
 *
 * @public
 * @returns {number}
 */
LocusZoom.DataLayer.prototype.getAbsoluteDataHeight = function() {
    var dataBCR = this.svg.group.node().getBoundingClientRect();
    return dataBCR.height;
};

/**
 * Whether transitions can be applied to this data layer
 * @returns {boolean}
 */
LocusZoom.DataLayer.prototype.canTransition = function() {
    if (!this.layout.transition) { return false; }
    return !(this.parent_plot.panel_boundaries.dragging || this.parent_plot.interaction.panel_id);
};

/**
 * Fetch the fully qualified ID to be associated with a specific visual element, based on the data to which that
 *   element is bound. In general this element ID will be unique, allowing it to be addressed directly via selectors.
 * @param {String|Object} element
 * @returns {String}
 */
LocusZoom.DataLayer.prototype.getElementId = function(element) {
    var element_id = 'element';
    if (typeof element == 'string') {
        element_id = element;
    } else if (typeof element == 'object') {
        var id_field = this.layout.id_field || 'id';
        if (typeof element[id_field] == 'undefined') {
            throw new Error('Unable to generate element ID');
        }
        element_id = element[id_field].toString().replace(/\W/g,'');
    }
    return (this.getBaseId() + '-' + element_id).replace(/([:.[\],])/g, '_');
};

/**
 * Fetch an ID that may bind a data element to a separate visual node for displaying status
 * Examples of this might be seperate visual nodes to show select/highlight statuses, or
 * even a common/shared node to show status across many elements in a set.
 * Abstract method. It should be overridden by data layers that implement seperate status
 * nodes specifically to the use case of the data layer type.
 * @param {String|Object} element
 * @returns {String|null}
 */
LocusZoom.DataLayer.prototype.getElementStatusNodeId = function(element) {
    return null;
};

/**
 * Returns a reference to the underlying data associated with a single visual element in the data layer, as
 *   referenced by the unique identifier for the element

 * @param {String} id The unique identifier for the element, as defined by `getElementId`
 * @returns {Object|null} The data bound to that element
 */
LocusZoom.DataLayer.prototype.getElementById = function(id) {
    var selector = d3.select('#' + id.replace(/([:.[\],])/g, '\\$1')); // escape special characters
    if (!selector.empty() && selector.data() && selector.data().length) {
        return selector.data()[0];
    } else {
        return null;
    }
};

/**
 * Basic method to apply arbitrary methods and properties to data elements.
 *   This is called on all data immediately after being fetched.
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.applyDataMethods = function() {
    this.data.forEach(function(d, i) {
        // Basic toHTML() method - return the stringified value in the id_field, if defined.
        this.data[i].toHTML = function() {
            var id_field = this.layout.id_field || 'id';
            var html = '';
            if (this.data[i][id_field]) { html = this.data[i][id_field].toString(); }
            return html;
        }.bind(this);
        // getDataLayer() method - return a reference to the data layer
        this.data[i].getDataLayer = function() {
            return this;
        }.bind(this);
        // deselect() method - shortcut method to deselect the element
        this.data[i].deselect = function() {
            var data_layer = this.getDataLayer();
            data_layer.unselectElement(this);
        };
    }.bind(this));
    this.applyCustomDataMethods();
    return this;
};

/**
 * Hook that allows custom datalayers to apply additional methods and properties to data elements as needed
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.applyCustomDataMethods = function() {
    return this;
};

/**
 * Initialize a data layer
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.initialize = function() {

    // Append a container group element to house the main data layer group element and the clip path
    this.svg.container = this.parent.svg.group.append('g')
        .attr('class', 'lz-data_layer-container')
        .attr('id', this.getBaseId() + '.data_layer_container');

    // Append clip path to the container element
    this.svg.clipRect = this.svg.container.append('clipPath')
        .attr('id', this.getBaseId() + '.clip')
        .append('rect');

    // Append svg group for rendering all data layer elements, clipped by the clip path
    this.svg.group = this.svg.container.append('g')
        .attr('id', this.getBaseId() + '.data_layer')
        .attr('clip-path', 'url(#' + this.getBaseId() + '.clip)');

    return this;

};

/**
 * Move a data layer up relative to others by z-index
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.moveUp = function() {
    if (this.parent.data_layer_ids_by_z_index[this.layout.z_index + 1]) {
        this.parent.data_layer_ids_by_z_index[this.layout.z_index] = this.parent.data_layer_ids_by_z_index[this.layout.z_index + 1];
        this.parent.data_layer_ids_by_z_index[this.layout.z_index + 1] = this.id;
        this.parent.resortDataLayers();
    }
    return this;
};

/**
 * Move a data layer down relative to others by z-index
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.moveDown = function() {
    if (this.parent.data_layer_ids_by_z_index[this.layout.z_index - 1]) {
        this.parent.data_layer_ids_by_z_index[this.layout.z_index] = this.parent.data_layer_ids_by_z_index[this.layout.z_index - 1];
        this.parent.data_layer_ids_by_z_index[this.layout.z_index - 1] = this.id;
        this.parent.resortDataLayers();
    }
    return this;
};

/**
 * Apply scaling functions to an element or parameter as needed, based on its layout and the element's data
 * If the layout parameter is already a primitive type, simply return the value as given
 * @param {Array|Number|String|Object} layout
 * @param {*} data The value to be used with the filter
 * @returns {*} The transformed value
 */
LocusZoom.DataLayer.prototype.resolveScalableParameter = function(layout, data) {
    var ret = null;
    if (Array.isArray(layout)) {
        var idx = 0;
        while (ret === null && idx < layout.length) {
            ret = this.resolveScalableParameter(layout[idx], data);
            idx++;
        }
    } else {
        switch (typeof layout) {
        case 'number':
        case 'string':
            ret = layout;
            break;
        case 'object':
            if (layout.scale_function) {
                if(layout.field) {
                    var f = new LocusZoom.Data.Field(layout.field);
                    ret = LocusZoom.ScaleFunctions.get(layout.scale_function, layout.parameters || {}, f.resolve(data));
                } else {
                    ret = LocusZoom.ScaleFunctions.get(layout.scale_function, layout.parameters || {}, data);
                }
            }
            break;
        }
    }
    return ret;
};


/**
 * Implementation hook for fetching the min and max values of available data. Used to determine axis range, if no other
 *   explicit axis settings override. Useful for data layers where the data extent depends on more than one field.
 *   (eg confidence intervals in a forest plot)
 * @param data
 * @param axis_config The configuration object for the specified axis.
 * @returns {Array} [min, max] without any padding applied
 * @private
 */
LocusZoom.DataLayer.prototype._getDataExtent = function(data, axis_config) {
    data = data || this.data;
    // By default this depends only on a single field.
    return d3.extent(data, function (d) {
        var f = new LocusZoom.Data.Field(axis_config.field);
        return +f.resolve(d);
    });
};

/**
 * Generate dimension extent function based on layout parameters
 * @param {('x'|'y')} dimension
 */
LocusZoom.DataLayer.prototype.getAxisExtent = function(dimension) {

    if (['x', 'y'].indexOf(dimension) === -1) {
        throw new Error('Invalid dimension identifier passed to LocusZoom.DataLayer.getAxisExtent()');
    }

    var axis_name = dimension + '_axis';
    var axis_layout = this.layout[axis_name];

    // If a floor AND a ceiling are explicitly defined then just return that extent and be done
    if (!isNaN(axis_layout.floor) && !isNaN(axis_layout.ceiling)) {
        return [+axis_layout.floor, +axis_layout.ceiling];
    }

    // If a field is defined for the axis and the data layer has data then generate the extent from the data set
    var data_extent = [];
    if (axis_layout.field && this.data) {
        if (!this.data.length) {
            // If data has been fetched (but no points in region), enforce the min_extent (with no buffers,
            //  because we don't need padding around an empty screen)
            data_extent = axis_layout.min_extent || [];
            return data_extent;
        } else {
            data_extent = this._getDataExtent(this.data, axis_layout);

            // Apply upper/lower buffers, if applicable
            var original_extent_span = data_extent[1] - data_extent[0];
            if (!isNaN(axis_layout.lower_buffer)) {
                data_extent[0] -= original_extent_span * axis_layout.lower_buffer;
            }
            if (!isNaN(axis_layout.upper_buffer)) {
                data_extent[1] += original_extent_span * axis_layout.upper_buffer;
            }

            if (typeof axis_layout.min_extent == 'object') {
                // The data should span at least the range specified by min_extent, an array with [low, high]
                var range_min = axis_layout.min_extent[0];
                var range_max = axis_layout.min_extent[1];
                if (!isNaN(range_min) && !isNaN(range_max)) {
                    data_extent[0] = Math.min(data_extent[0], range_min);
                }
                if (!isNaN(range_max)) {
                    data_extent[1] = Math.max(data_extent[1], range_max);
                }
            }
            // If specified, floor and ceiling will override the actual data range
            return [
                isNaN(axis_layout.floor) ? data_extent[0] : axis_layout.floor,
                isNaN(axis_layout.ceiling) ? data_extent[1] : axis_layout.ceiling
            ];
        }
    }

    // If this is for the x axis and no extent could be generated yet but state has a defined start and end
    // then default to using the state-defined region as the extent
    if (dimension === 'x' && !isNaN(this.state.start) && !isNaN(this.state.end)) {
        return [this.state.start, this.state.end];
    }

    // No conditions met for generating a valid extent, return an empty array
    return [];

};

/**
 * Allow this data layer to tell the panel what axis ticks it thinks it will require. The panel may choose whether
 *   to use some, all, or none of these when rendering, either alone or in conjunction with other data layers.
 *
 *   This method is a stub and should be overridden in data layers that need to specify custom behavior.
 *
 * @param {('x'|'y1'|'y2')} dimension
 * @param {Object} [config] Additional parameters for the panel to specify how it wants ticks to be drawn. The names
 *   and meanings of these parameters may vary between different data layers.
 * @returns {Object[]}
 *   An array of objects: each object must have an 'x' attribute to position the tick.
 *   Other supported object keys:
 *     * text: string to render for a given tick
 *     * style: d3-compatible CSS style object
 *     * transform: SVG transform attribute string
 *     * color: string or LocusZoom scalable parameter object
 */
LocusZoom.DataLayer.prototype.getTicks = function (dimension, config) {
    if (['x', 'y1', 'y2'].indexOf(dimension) === -1) {
        throw new Error('Invalid dimension identifier at layer level' + dimension);
    }
    return [];
};

/**
 * Generate a tool tip for a given element
 * @param {String|Object} d The element associated with the tooltip
 * @param {String} [id] An identifier to the tooltip
 */
LocusZoom.DataLayer.prototype.createTooltip = function(d, id) {
    if (typeof this.layout.tooltip != 'object') {
        throw new Error('DataLayer [' + this.id + '] layout does not define a tooltip');
    }
    if (typeof id == 'undefined') { id = this.getElementId(d); }
    if (this.tooltips[id]) {
        this.positionTooltip(id);
        return;
    }
    this.tooltips[id] = {
        data: d,
        arrow: null,
        selector: d3.select(this.parent_plot.svg.node().parentNode).append('div')
            .attr('class', 'lz-data_layer-tooltip')
            .attr('id', id + '-tooltip')
    };
    this.updateTooltip(d);
    return this;
};

/**
 * Update a tool tip (generate its inner HTML)
 * @param {String|Object} d The element associated with the tooltip
 * @param {String} [id] An identifier to the tooltip
 */
LocusZoom.DataLayer.prototype.updateTooltip = function(d, id) {
    if (typeof id == 'undefined') { id = this.getElementId(d); }
    // Empty the tooltip of all HTML (including its arrow!)
    this.tooltips[id].selector.html('');
    this.tooltips[id].arrow = null;
    // Set the new HTML
    if (this.layout.tooltip.html) {
        this.tooltips[id].selector.html(LocusZoom.parseFields(d, this.layout.tooltip.html));
    }
    // If the layout allows tool tips on this data layer to be closable then add the close button
    // and add padding to the tooltip to accommodate it
    if (this.layout.tooltip.closable) {
        this.tooltips[id].selector.insert('button', ':first-child')
            .attr('class', 'lz-tooltip-close-button')
            .attr('title', 'Close')
            .text('×')
            .on('click', function() {
                this.destroyTooltip(id);
            }.bind(this));
    }
    // Apply data directly to the tool tip for easier retrieval by custom UI elements inside the tool tip
    this.tooltips[id].selector.data([d]);
    // Reposition and draw a new arrow
    this.positionTooltip(id);
    return this;
};

/**
 * Destroy tool tip - remove the tool tip element from the DOM and delete the tool tip's record on the data layer
 * @param {String|Object} d The element associated with the tooltip
 * @param {String} [id] An identifier to the tooltip
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.destroyTooltip = function(d, id) {
    if (typeof d == 'string') {
        id = d;
    } else if (typeof id == 'undefined') {
        id = this.getElementId(d);
    }
    if (this.tooltips[id]) {
        if (typeof this.tooltips[id].selector == 'object') {
            this.tooltips[id].selector.remove();
        }
        delete this.tooltips[id];
    }
    return this;
};

/**
 * Loop through and destroy all tool tips on this data layer
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.destroyAllTooltips = function() {
    for (var id in this.tooltips) {
        this.destroyTooltip(id);
    }
    return this;
};

//
/**
 * Position tool tip - naïve function to place a tool tip to the lower right of the current mouse element
 *   Most data layers reimplement this method to position tool tips specifically for the data they display
 * @param {String} id The identifier of the tooltip to position
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.positionTooltip = function(id) {
    if (typeof id != 'string') {
        throw new Error('Unable to position tooltip: id is not a string');
    }
    // Position the div itself
    this.tooltips[id].selector
        .style('left', (d3.event.pageX) + 'px')
        .style('top', (d3.event.pageY) + 'px');
    // Create / update position on arrow connecting tooltip to data
    if (!this.tooltips[id].arrow) {
        this.tooltips[id].arrow = this.tooltips[id].selector.append('div')
            .style('position', 'absolute')
            .attr('class', 'lz-data_layer-tooltip-arrow_top_left');
    }
    this.tooltips[id].arrow
        .style('left', '-1px')
        .style('top', '-1px');
    return this;
};

/**
 * Loop through and position all tool tips on this data layer
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.positionAllTooltips = function() {
    for (var id in this.tooltips) {
        this.positionTooltip(id);
    }
    return this;
};

/**
 * Show or hide a tool tip by ID depending on directives in the layout and state values relative to the ID
 * @param {String|Object} element The element associated with the tooltip
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.showOrHideTooltip = function(element) {

    if (typeof this.layout.tooltip != 'object') { return; }
    var id = this.getElementId(element);

    var resolveStatus = function(statuses, directive, operator) {
        var status = null;
        if (typeof statuses != 'object' || statuses === null) { return null; }
        if (Array.isArray(directive)) {
            if (typeof operator == 'undefined') { operator = 'and'; }
            if (directive.length === 1) {
                status = statuses[directive[0]];
            } else {
                status = directive.reduce(function(previousValue, currentValue) {
                    if (operator === 'and') {
                        return statuses[previousValue] && statuses[currentValue];
                    } else if (operator === 'or') {
                        return statuses[previousValue] || statuses[currentValue];
                    }
                    return null;
                });
            }
        } else if (typeof directive == 'object') {
            var sub_status;
            for (var sub_operator in directive) {
                sub_status = resolveStatus(statuses, directive[sub_operator], sub_operator);
                if (status === null) {
                    status = sub_status;
                } else if (operator === 'and') {
                    status = status && sub_status;
                } else if (operator === 'or') {
                    status = status || sub_status;
                }
            }
        }
        return status;
    };

    var show_directive = {};
    if (typeof this.layout.tooltip.show == 'string') {
        show_directive = { and: [ this.layout.tooltip.show ] };
    } else if (typeof this.layout.tooltip.show == 'object') {
        show_directive = this.layout.tooltip.show;
    }

    var hide_directive = {};
    if (typeof this.layout.tooltip.hide == 'string') {
        hide_directive = { and: [ this.layout.tooltip.hide ] };
    } else if (typeof this.layout.tooltip.hide == 'object') {
        hide_directive = this.layout.tooltip.hide;
    }

    var statuses = {};
    LocusZoom.DataLayer.Statuses.adjectives.forEach(function(status) {
        var antistatus = 'un' + status;
        statuses[status] = this.state[this.state_id][status].indexOf(id) !== -1;
        statuses[antistatus] = !statuses[status];
    }.bind(this));

    var show_resolved = resolveStatus(statuses, show_directive);
    var hide_resolved = resolveStatus(statuses, hide_directive);

    // Only show tooltip if the resolved logic explicitly shows and explicitly not hides the tool tip
    // Otherwise ensure tooltip does not exist
    if (show_resolved && !hide_resolved) {
        this.createTooltip(element);
    } else {
        this.destroyTooltip(element);
    }

    return this;

};

/**
 * Find the elements (or indices) that match any of a set of provided filters
 * @protected
 * @param {Array[]} filters A list of filter entries: [field, value] (for equivalence testing) or
 *   [field, operator, value] for other operators
 * @param {('indexes'|'elements')} [return_type='indexes'] Specify whether to return either the indices of the matching
 *   elements, or references to the elements themselves
 * @returns {Array}
 */
LocusZoom.DataLayer.prototype.filter = function(filters, return_type) {
    if (typeof return_type == 'undefined' || ['indexes','elements'].indexOf(return_type) === -1) {
        return_type = 'indexes';
    }
    if (!Array.isArray(filters)) { return []; }
    var test = function(element, filter) {
        var operators = {
            '=': function(a,b) { return a === b; },
            // eslint-disable-next-line eqeqeq
            '!=': function(a,b) { return a != b; }, // For absence of a value, deliberately allow weak comparisons (eg undefined/null)
            '<': function(a,b) { return a < b; },
            '<=': function(a,b) { return a <= b; },
            '>': function(a,b) { return a > b; },
            '>=': function(a,b) { return a >= b; },
            '%': function(a,b) { return a % b; }
        };
        if (!Array.isArray(filter)) { return false; }
        if (filter.length === 2) {
            return element[filter[0]] === filter[1];
        } else if (filter.length === 3 && operators[filter[1]]) {
            return operators[filter[1]](element[filter[0]], filter[2]);
        } else {
            return false;
        }
    };
    var matches = [];
    this.data.forEach(function(element, idx) {
        var match = true;
        filters.forEach(function(filter) {
            if (!test(element, filter)) { match = false; }
        });
        if (match) { matches.push(return_type === 'indexes' ? idx : element); }
    });
    return matches;
};

/**
 * @param filters
 * @returns {Array}
 */
LocusZoom.DataLayer.prototype.filterIndexes = function(filters) { return this.filter(filters, 'indexes'); };
/**
 * @param filters
 * @returns {Array}
 */
LocusZoom.DataLayer.prototype.filterElements = function(filters) { return this.filter(filters, 'elements'); };

LocusZoom.DataLayer.Statuses.verbs.forEach(function(verb, idx) {
    var adjective = LocusZoom.DataLayer.Statuses.adjectives[idx];
    var antiverb = 'un' + verb;
    // Set/unset a single element's status
    // TODO: Improve documentation for dynamically generated methods/properties
    LocusZoom.DataLayer.prototype[verb + 'Element'] = function(element, exclusive) {
        if (typeof exclusive == 'undefined') { exclusive = false; } else { exclusive = !!exclusive; }
        this.setElementStatus(adjective, element, true, exclusive);
        return this;
    };
    LocusZoom.DataLayer.prototype[antiverb + 'Element'] = function(element, exclusive) {
        if (typeof exclusive == 'undefined') { exclusive = false; } else { exclusive = !!exclusive; }
        this.setElementStatus(adjective, element, false, exclusive);
        return this;
    };
    // Set/unset status for arbitrarily many elements given a set of filters
    LocusZoom.DataLayer.prototype[verb + 'ElementsByFilters'] = function(filters, exclusive) {
        if (typeof exclusive == 'undefined') { exclusive = false; } else { exclusive = !!exclusive; }
        return this.setElementStatusByFilters(adjective, true, filters, exclusive);
    };
    LocusZoom.DataLayer.prototype[antiverb + 'ElementsByFilters'] = function(filters, exclusive) {
        if (typeof exclusive == 'undefined') { exclusive = false; } else { exclusive = !!exclusive; }
        return this.setElementStatusByFilters(adjective, false, filters, exclusive);
    };
    // Set/unset status for all elements
    LocusZoom.DataLayer.prototype[verb + 'AllElements'] = function() {
        this.setAllElementStatus(adjective, true);
        return this;
    };
    LocusZoom.DataLayer.prototype[antiverb + 'AllElements'] = function() {
        this.setAllElementStatus(adjective, false);
        return this;
    };
});

/**
 * Toggle a status (e.g. highlighted, selected, identified) on an element
 * @param {String} status The name of a recognized status to be added/removed on an appropriate element
 * @param {String|Object} element The data bound to the element of interest
 * @param {Boolean} toggle True to add the status (and associated CSS styles); false to remove it
 * @param {Boolean} exclusive Whether to only allow a state for a single element at a time
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.setElementStatus = function(status, element, toggle, exclusive) {
    // Sanity checks
    if (typeof status == 'undefined' || LocusZoom.DataLayer.Statuses.adjectives.indexOf(status) === -1) {
        throw new Error('Invalid status passed to DataLayer.setElementStatus()');
    }
    if (typeof element == 'undefined') {
        throw new Error('Invalid element passed to DataLayer.setElementStatus()');
    }
    if (typeof toggle == 'undefined') {
        toggle = true;
    }

    // Get an ID for the element or return having changed nothing
    try {
        var element_id = this.getElementId(element);
    } catch (get_element_id_error) {
        return this;
    }

    // Enforce exclusivity (force all elements to have the opposite of toggle first)
    if (exclusive) {
        this.setAllElementStatus(status, !toggle);
    }

    // Set/unset the proper status class on the appropriate DOM element(s)
    d3.select('#' + element_id).classed('lz-data_layer-' + this.layout.type + '-' + status, toggle);
    var element_status_node_id = this.getElementStatusNodeId(element);
    if (element_status_node_id !== null) {
        d3.select('#' + element_status_node_id).classed('lz-data_layer-' + this.layout.type + '-statusnode-' + status, toggle);
    }

    // Track element ID in the proper status state array
    var element_status_idx = this.state[this.state_id][status].indexOf(element_id);
    if (toggle && element_status_idx === -1) {
        this.state[this.state_id][status].push(element_id);
    }
    if (!toggle && element_status_idx !== -1) {
        this.state[this.state_id][status].splice(element_status_idx, 1);
    }

    // Trigger tool tip show/hide logic
    this.showOrHideTooltip(element);

    // Trigger layout changed event hook
    this.parent.emit('layout_changed', true);
    if (status === 'selected') {
        // Notify parents that a given element has been interacted with. For now, we will only notify on
        //   "selected" type events, which is (usually) a toggle-able state. If elements are exclusive, two selection
        //   events will be sent in short order as the previously selected element has to be de-selected first
        this.parent.emit('element_selection', { element: element, active: toggle }, true);
    }
    return this;

};

/**
 * Toggle a status on elements in the data layer based on a set of filters
 * @param {String} status
 * @param {Boolean} toggle
 * @param {Array} filters
 * @param {Boolean} exclusive
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.setElementStatusByFilters = function(status, toggle, filters, exclusive) {

    // Sanity check
    if (typeof status == 'undefined' || LocusZoom.DataLayer.Statuses.adjectives.indexOf(status) === -1) {
        throw new Error('Invalid status passed to DataLayer.setElementStatusByFilters()');
    }
    if (typeof this.state[this.state_id][status] == 'undefined') { return this; }
    if (typeof toggle == 'undefined') { toggle = true; } else { toggle = !!toggle; }
    if (typeof exclusive == 'undefined') { exclusive = false; } else { exclusive = !!exclusive; }
    if (!Array.isArray(filters)) { filters = []; }

    // Enforce exclusivity (force all elements to have the opposite of toggle first)
    if (exclusive) {
        this.setAllElementStatus(status, !toggle);
    }

    // Apply statuses
    this.filterElements(filters).forEach(function(element) {
        this.setElementStatus(status, element, toggle);
    }.bind(this));

    return this;
};

/**
 * Toggle a status on all elements in the data layer
 * @param {String} status
 * @param {Boolean} toggle
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.setAllElementStatus = function(status, toggle) {

    // Sanity check
    if (typeof status == 'undefined' || LocusZoom.DataLayer.Statuses.adjectives.indexOf(status) === -1) {
        throw new Error('Invalid status passed to DataLayer.setAllElementStatus()');
    }
    if (typeof this.state[this.state_id][status] == 'undefined') { return this; }
    if (typeof toggle == 'undefined') { toggle = true; }

    // Apply statuses
    if (toggle) {
        this.data.forEach(function(element) {
            this.setElementStatus(status, element, true);
        }.bind(this));
    } else {
        var status_ids = this.state[this.state_id][status].slice();
        status_ids.forEach(function(id) {
            var element = this.getElementById(id);
            if (typeof element == 'object' && element !== null) {
                this.setElementStatus(status, element, false);
            }
        }.bind(this));
        this.state[this.state_id][status] = [];
    }

    // Update global status flag
    this.global_statuses[status] = toggle;

    return this;
};

/**
 * Apply all layout-defined behaviors (DOM event handlers) to a selection of elements
 * @param {d3.selection} selection
 */
LocusZoom.DataLayer.prototype.applyBehaviors = function(selection) {
    if (typeof this.layout.behaviors != 'object') { return; }
    Object.keys(this.layout.behaviors).forEach(function(directive) {
        var event_match = /(click|mouseover|mouseout)/.exec(directive);
        if (!event_match) { return; }
        selection.on(event_match[0] + '.' + directive, this.executeBehaviors(directive, this.layout.behaviors[directive]));
    }.bind(this));
};

/**
 * Generate a function that executes an arbitrary list of behaviors on an element during an event
 * @param {String} directive The name of the event, as described in layout.behaviors for this datalayer
 * @param {Object} behaviors An object describing the behavior to attach to this single element
 * @param {string} behaviors.action The name of the action that would trigger this behavior (eg click, mouseover, etc)
 * @param {string} behaviors.status What status to apply to the element when this behavior is triggered (highlighted,
 *  selected, etc)
 * @param {string} [behaviors.exclusive] Whether triggering the event for this element should unset the relevant status
 *   for all other elements. Useful for, eg, click events that exclusively highlight one thing.
 * @returns {function(this:LocusZoom.DataLayer)} Return a function that handles the event in context with the behavior
 *   and the element- can be attached as an event listener
 */
LocusZoom.DataLayer.prototype.executeBehaviors = function(directive, behaviors) {

    // Determine the required state of control and shift keys during the event
    var requiredKeyStates = {
        'ctrl': (directive.indexOf('ctrl') !== -1),
        'shift': (directive.indexOf('shift') !== -1)
    };

    return function(element) {

        // Do nothing if the required control and shift key presses (or lack thereof) doesn't match the event
        if (requiredKeyStates.ctrl !== !!d3.event.ctrlKey || requiredKeyStates.shift !== !!d3.event.shiftKey) { return; }

        // Loop through behaviors making each one go in succession
        behaviors.forEach(function(behavior) {

            // Route first by the action, if defined
            if (typeof behavior != 'object' || behavior === null) { return; }

            switch (behavior.action) {

            // Set a status (set to true regardless of current status, optionally with exclusivity)
            case 'set':
                this.setElementStatus(behavior.status, element, true, behavior.exclusive);
                break;

            // Unset a status (set to false regardless of current status, optionally with exclusivity)
            case 'unset':
                this.setElementStatus(behavior.status, element, false, behavior.exclusive);
                break;

            // Toggle a status
            case 'toggle':
                var current_status_boolean = (this.state[this.state_id][behavior.status].indexOf(this.getElementId(element)) !== -1);
                var exclusive = behavior.exclusive && !current_status_boolean;
                this.setElementStatus(behavior.status, element, !current_status_boolean, exclusive);
                break;

            // Link to a dynamic URL
            case 'link':
                if (typeof behavior.href == 'string') {
                    var url = LocusZoom.parseFields(element, behavior.href);
                    if (typeof behavior.target == 'string') {
                        window.open(url, behavior.target);
                    } else {
                        window.location.href = url;
                    }
                }
                break;

            // Action not defined, just return
            default:
                break;

            }

            return;

        }.bind(this));

    }.bind(this);

};

/**
 * Get an object with the x and y coordinates of the panel's origin in terms of the entire page
 *   Necessary for positioning any HTML elements over the panel
 * @returns {{x: Number, y: Number}}
 */
LocusZoom.DataLayer.prototype.getPageOrigin = function() {
    var panel_origin = this.parent.getPageOrigin();
    return {
        x: panel_origin.x + this.parent.layout.margin.left,
        y: panel_origin.y + this.parent.layout.margin.top
    };
};

/**
 * Get a data layer's current underlying data in a standard format (e.g. JSON or CSV)
 * @param {('csv'|'tsv'|'json')} format How to export the data
 * @returns {*}
 */
LocusZoom.DataLayer.prototype.exportData = function(format) {
    var default_format = 'json';
    format = format || default_format;
    format = (typeof format == 'string' ? format.toLowerCase() : default_format);
    if (['json','csv','tsv'].indexOf(format) === -1) { format = default_format; }
    var ret;
    switch (format) {
    case 'json':
        try {
            ret = JSON.stringify(this.data);
        } catch (e) {
            ret = null;
            console.warn('Unable to export JSON data from data layer: ' + this.getBaseId());
            console.error(e);
        }
        break;
    case 'tsv':
    case 'csv':
        try {
            var jsonified = JSON.parse(JSON.stringify(this.data));
            if (typeof jsonified != 'object') {
                ret = jsonified.toString();
            } else if (!Array.isArray(jsonified)) {
                ret = 'Object';
            } else {
                var delimiter = (format === 'tsv') ? '\t' : ',';
                var header = this.layout.fields.map(function(header) {
                    return JSON.stringify(header);
                }).join(delimiter) + '\n';
                ret = header + jsonified.map(function(record) {
                    return this.layout.fields.map(function(field) {
                        if (typeof record[field] == 'undefined') {
                            return JSON.stringify(null);
                        } else if (typeof record[field] == 'object' && record[field] !== null) {
                            return Array.isArray(record[field]) ? '"[Array(' + record[field].length + ')]"' : '"[Object]"';
                        } else {
                            return JSON.stringify(record[field]);
                        }
                    }).join(delimiter);
                }.bind(this)).join('\n');
            }
        } catch (e) {
            ret = null;
            console.error('Unable to export CSV data from data layer: ' + this.getBaseId() + ';', e);
        }
        break;
    }
    return ret;
};

/**
 * Position the datalayer and all tooltips
 * @returns {LocusZoom.DataLayer}
 */
LocusZoom.DataLayer.prototype.draw = function() {
    this.svg.container.attr('transform', 'translate(' + this.parent.layout.cliparea.origin.x +  ',' + this.parent.layout.cliparea.origin.y + ')');
    this.svg.clipRect
        .attr('width', this.parent.layout.cliparea.width)
        .attr('height', this.parent.layout.cliparea.height);
    this.positionAllTooltips();
    return this;
};


/**
 * Re-Map a data layer to reflect changes in the state of a plot (such as viewing region/ chromosome range)
 * @return {Promise}
 */
LocusZoom.DataLayer.prototype.reMap = function() {

    this.destroyAllTooltips(); // hack - only non-visible tooltips should be destroyed
    // and then recreated if returning to visibility

    // Fetch new data. Datalayers are only given access to the final consolidated data from the chain (not headers or raw payloads)
    var promise = this.parent_plot.lzd.getData(this.state, this.layout.fields);
    promise.then(function(new_data) {
        this.data = new_data.body;
        this.applyDataMethods();
        this.initialized = true;
    }.bind(this));

    return promise;
};


/**
 * The central registry of known data layer definitions (which may be stored in separate files due to length)
 * @namespace
 */
LocusZoom.DataLayers = (function() {
    var obj = {};
    var datalayers = {};
    /**
     * @name LocusZoom.DataLayers.get
     * @param {String} name The name of the datalayer
     * @param {Object} layout The configuration object for this data layer
     * @param {LocusZoom.DataLayer|LocusZoom.Panel} parent Where this layout is used
     * @returns {LocusZoom.DataLayer}
     */
    obj.get = function(name, layout, parent) {
        if (!name) {
            return null;
        } else if (datalayers[name]) {
            if (typeof layout != 'object') {
                throw new Error('invalid layout argument for data layer [' + name + ']');
            } else {
                return new datalayers[name](layout, parent);
            }
        } else {
            throw new Error('data layer [' + name + '] not found');
        }
    };

    /**
     * @name LocusZoom.DataLayers.set
     * @protected
     * @param {String} name
     * @param {Function} datalayer Constructor for the datalayer
     */
    obj.set = function(name, datalayer) {
        if (datalayer) {
            if (typeof datalayer != 'function') {
                throw new Error('unable to set data layer [' + name + '], argument provided is not a function');
            } else {
                datalayers[name] = datalayer;
                datalayers[name].prototype = new LocusZoom.DataLayer();
            }
        } else {
            delete datalayers[name];
        }
    };

    /**
     * Add a new type of datalayer to the registry of known layer types
     * @name LocusZoom.DataLayers.add
     * @param {String} name The name of the data layer to register
     * @param {Function} datalayer
     */
    obj.add = function(name, datalayer) {
        if (datalayers[name]) {
            throw new Error('data layer already exists with name: ' + name);
        } else {
            obj.set(name, datalayer);
        }
    };

    /**
     * Register a new datalayer that inherits and extends basic behaviors from a known datalayer
     * @param {String} parent_name The name of the parent data layer whose behavior is to be extended
     * @param {String} name The name of the new datalayer to register
     * @param {Object} [overrides] Object of properties and methods to combine with the prototype of the parent datalayer
     * @returns {Function} The constructor for the new child class
     */
    obj.extend = function(parent_name, name, overrides) {
        // TODO: Consider exposing additional constructor argument, if there is a use case for very granular extension
        overrides = overrides || {};

        var parent = datalayers[parent_name];
        if (!parent) {
            throw new Error('Attempted to subclass an unknown or unregistered datalayer type');
        }
        if (typeof overrides !== 'object') {
            throw new Error('Must specify an object of properties and methods');
        }
        var child = LocusZoom.subclass(parent, overrides);
        // Bypass .set() because we want a layer of inheritance below `DataLayer`
        datalayers[name] = child;
        return child;
    };

    /**
     * List the names of all known datalayers
     * @name LocusZoom.DataLayers.list
     * @returns {String[]}
     */
    obj.list = function() {
        return Object.keys(datalayers);
    };

    return obj;
})();
